import os
import re
import subprocess
import uuid
from collections import defaultdict

import bam_utils

ANNOTATION_TO_INDICES = {
        'qseqid': 0,
        'sseqid': 1,
        'pident': 2,
        'qlen': 3,
        'qstart': 6,
        'qend': 7,
        'sstart': 8,
        'send': 9,
        'evalue': 10,
        'bitscore': 11
        }

def parse_blat_output(output_fp):
    """Returns list of dicts representing each line in output"""
    output_dicts = []
    f = open(output_fp)
    for line in f:
        d = {}
        pieces = line.strip().split('\t')
        if pieces:
            for field, index in ANNOTATION_TO_INDICES.items():
                d[field] = pieces[index]

            output_dicts.append(d)
    f.close()

    return output_dicts

def execute_blat(query_fp, database, out='blast8', output_fp='temp.out'):
    tool_args = ['blat', database, query_fp,
            f'-out={out}',
            output_fp]

    subprocess.check_output(tool_args).decode('utf-8')

def is_in_range(chrom, read_start, read_end, d):
    norm_chrom = re.sub(r'^chr(.*)$', r'\1', chrom)
    chrom = re.sub(r'^chr(.*)$', r'\1', d['sseqid'])
    start_pos = int(d['sstart'])
    end_pos = int(d['send'])

    start_pos_in_range = start_pos <= read_end + 1 and start_pos >= read_start - 1
    end_pos_in_range = end_pos <= read_end + 1 and end_pos >= read_start - 1
    if chrom == norm_chrom and start_pos_in_range and end_pos_in_range:
        return True
    return False

def is_positive_rna_count(chrom, read_start, read_end, blat_result_dicts, percent_threshold=.95):
    norm_chrom = re.sub(r'^chr(.*)$', r'\1', chrom)

    # sort blast result dicts by score
    blat_result_dicts = sorted(blat_result_dicts, key=lambda x: float(x['bitscore']), reverse=True)

    # if no passing return false
    if len(blat_result_dicts) == 0:
        return False

    # check first entry to see if it is in right position
    d = blat_result_dicts[0]
    if not is_in_range(chrom, read_start, read_end, d):
        return False
    score = float(d['bitscore'])

    # if only dict return true
    if len(blat_result_dicts) == 1:
        return True

    # check second entry to see if it meets threshold
    d = blat_result_dicts[1]
    if float(d['bitscore']) > percent_threshold * score:
        return False

    return True

class BlatAnnotator(object):
    def __init__(self, annotations, database, rna_editing_percent_threshold=.95):
        self.annotations = annotations
        self.database = database

        self.rna_editing_percent_threshold = rna_editing_percent_threshold

        self.reads_to_data = {}
        self.position_to_reference_base = {}

    def prepare_input_files(self, input_bam_fp, output_fasta_fp, position_tups):
        """prepare input files that BlastAnnotator needs if reading from bam and position file
    
        position_tups - [(chrom, pos), ...]"""
        # index the bam in case it isn't already
        bam_utils.index_bam(input_bam_fp)

        # create a positions file that will work with samtools in case input doesn't
        u_id = str(uuid.uuid4())
        temp_positions_fp = f'temp.positions.{u_id}.bed'
        out_f = open(temp_positions_fp, 'w')
        for chrom, pos in position_tups:
            out_f.write(f'{chrom}\t{pos}\t{pos}\n')
        out_f.close()

        self.reads_to_data = bam_utils.write_position_fasta(input_bam_fp,
                temp_positions_fp, output_fasta_fp)

        os.remove(temp_positions_fp)

    def blat_fasta(self, input_fasta):
        """Blat the given fasta and collect results for each sequence in input fasta"""
        u_id = str(uuid.uuid4())
        temp_output_fp = f'temp.{u_id}.out'
        execute_blat(input_fasta, self.database, output_fp=temp_output_fp)
        output_dicts = parse_blat_output(temp_output_fp)

        # remove temp output
        os.remove(temp_output_fp)

        sequence_to_results = defaultdict(list)
        for d in output_dicts:
            if 'qseqid' in d:
                sequence_to_results[d['qseqid']].append(d)

        return sequence_to_results

    def get_rna_editing_annotations(self, position_to_read_results):
        """Returns dict {position: %passing}"""
        position_to_percent_passing = {}
        for (chrom, pos), read_to_result_dicts in position_to_read_results.items():
            count, total = 0, 0
            for read, result_dicts in read_to_result_dicts.items():
                read_data = self.reads_to_data[f'{chrom}:{pos}|{read}']
                
                reference_base = self.position_to_reference_base[(chrom, str(pos))]
                read_start, read_end = bam_utils.get_covering_reference_coords(int(read_data['start']),
                        read_data['cigar'], read_data['sequence'])
                read_base = bam_utils.get_base_by_position(int(read_data['start']), int(pos),
                        read_data['cigar'], read_data['sequence'])

                if reference_base is not None and read_base is not None:
                    if reference_base.lower() != read_base.lower():
                        total += 1
                        if is_positive_rna_count(chrom, read_start, read_end, result_dicts,
                                percent_threshold=self.rna_editing_percent_threshold):
                            count += 1

            position_to_percent_passing[(chrom, pos)] = count / max(1, len(read_to_result_dicts))

        return position_to_percent_passing

    def get_rna_editing_blat_annotations(self, input_fasta,
            chrom_regex=r'^(.*):.*$', pos_regex=r'^.*:(.*)\|.*$', read_regex=r'^.*\|(.*)$'):
        """Collect by value in input fasta sequence identifier. Identifier is seperated by |

        i.e. if a sequence id is chr1:12345|read1, then the returned dictionary will
        look something like - {'chr1:12345': {read1: [{blastn parsed result}, ...], ...}, ...}
        """
        sequence_to_results = self.blat_fasta(input_fasta)

        position_to_read_results = {}
        for sequence_id, result_dict in sequence_to_results.items():
            chrom = re.sub(chrom_regex, r'\1', sequence_id)
            pos = int(re.sub(pos_regex, r'\1', sequence_id))
            read = re.sub(read_regex, r'\1', sequence_id)
            
            pos_tup = (chrom, pos)

            if pos_tup not in position_to_read_results:
                position_to_read_results[pos_tup] = {}
            position_to_read_results[pos_tup][read] = result_dict

        position_to_percent_passing = self.get_rna_editing_annotations(position_to_read_results)

        return position_to_percent_passing

    def get_blat_annotations_for_bam(self, input_bam_fp, position_tups, reference_bases=None):
        """Get annotations for the given positions based on reads in the given bam.

        input_bam_fp - filepath to bam with reads covering positions in the positions file
        position_tups - positions to recieve annotations. format - [(chrom, pos), ...]
        reference_bases - 

        if position is not present in bam, then it will not be returned in annotations

        Returns: position_to_annotations, headers
            {(chrom, pos): [annotation1, annotation2, annotation3, ...]}, [header1, 
                    header2, header3, ...]
        """
        u_id = str(uuid.uuid4())
        temp_fasta_fp = f'temp.query.{u_id}.fa'
        self.prepare_input_files(input_bam_fp, temp_fasta_fp, position_tups) 

        annotations_dict = defaultdict(list)
        headers = []
        if 'rna_editing' in self.annotations:
            if reference_bases is None:
                raise ValueError('reference bases must be present if doing rna editing annotations')

            self.position_to_reference_base = {(chrom, pos):base
                    for (chrom, pos), base in zip(position_tups, reference_bases)}

            position_to_percent_passing = self.get_rna_editing_blat_annotations(temp_fasta_fp)
            # add positions that were missing for whatever reason
            for (chrom, pos) in position_tups:
                if (chrom, int(pos)) not in position_to_percent_passing:
                    position_to_percent_passing[(chrom, int(pos))] = '.'

            headers += ['BLAT_RNA_EDITING_%_PASSING']
            for (chrom, pos), value in position_to_percent_passing.items():
                annotations_dict[(chrom, str(pos))].append(value)

        os.remove(temp_fasta_fp)

        return annotations_dict, headers
