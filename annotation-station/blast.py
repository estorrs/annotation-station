import os
import re
import subprocess
from collections import defaultdict

import bam_utils

ANNOTATION_TO_OUT_FIELDS = {
        'rna_editing': ['qseqid', 'sseqid', 'sstart', 'send', 'qstart', 'qend', 'qlen',
                'evalue', 'bitscore', 'pident'],
        'avg_top_hits': ['qseqid', 'sseqid', 'evalue', 'bitscore', 'pident']
        }

def get_required_out_fields(annotations):
    """return all needed output fields for the given annotations"""
    out_fields = set()
    for annotation in annotations:
        out_fields.update(ANNOTATION_TO_OUT_FIELDS[annotation])

    return list(out_fields)

def parse_blast_output(output, out_fields):
    """Returns list of dicts representing each line in output"""
    output_dicts = []
    lines = output.split('\n')
    for line in lines:
        d = {}
        pieces = line.split('\t')
        for field, val in zip(out_fields, pieces):
            d[field] = val

        output_dicts.append(d)

    return output_dicts

def execute_blastn(query_fp, database,
        out_fields=['qseqid', 'sseqid', 'sstart', 'send', 'qstart', 'qend', 'evalue', 'bitscore', 'pident'],
        max_target_seqs=5, max_hsps=5):
    outfmt_str = ' '.join(out_fields)
    tool_args = ['blastn',
            '-query', query_fp,
            '-db', database,
            '-task', 'blastn',
            '-max_target_seqs', str(max_target_seqs),
            '-max_hsps', str(max_hsps),
            '-outfmt', f'6 {outfmt_str}']

    result = subprocess.check_output(tool_args).decode('utf-8')

    return result

def is_in_range(chrom, pos, d):
    norm_chrom = re.sub(r'^chr(.*)$', r'\1', chrom)
    chrom = re.sub(r'^chr(.*)$', r'\1', d['sseqid'])
    start_pos = int(d['sstart'])
    end_pos = int(d['send'])
    if chrom == norm_chrom and start_pos <= pos and end_pos >= pos:
        return True
    return False

def is_positive_rna_count(chrom, pos, blast_result_dicts,
        identity_threshold=.95, coverage_threshold=.9):
    norm_chrom = re.sub(r'^chr(.*)$', r'\1', chrom)

    # sort blast result dicts by score
    blast_result_dicts = sorted(blast_result_dicts, key=lambda x: float(x['bitscore']), reverse=True)

    # throw out less than 90% coverage
    blast_result_dicts = [d for d in blast_result_dicts
            if (int(d['qend']) - int(d['qstart'])) / int(d['qlen']) >= coverage_threshold]

    # if no passing return false
    if len(blast_result_dicts) == 0:
        return False

    # check first entry to see if it is in right position and meets threshold
    d = blast_result_dicts[0]
    if not is_in_range(chrom, pos, d) or float(d['pident']) / 100 < identity_threshold:
        return False

    # if only dict return true
    if len(blast_result_dicts) == 1:
        return True

    # check second entry to see if it meets threshold
    d = blast_result_dicts[1]
    if float(d['pident']) / 100 > identity_threshold:
        return False

    return True

def prepare_input_files(input_bam_fp, output_fasta_fp, position_tups):
    """prepare input files that BlastAnnotator needs if reading from bam and position file
    
    position_tups - [(chrom, pos), ...]"""
    # index the bam in case it isn't already
    bam_utils.index_bam(input_bam_fp)

    # create a positions file that will work with samtools in case input doesn't
    temp_positions_fp = 'temp.positions.bed'
    out_f = open(temp_positions_fp, 'w')
    for chrom, pos in position_tups:
        out_f.write(f'{chrom}\t{pos}\t{pos}\n')
    out_f.close()

    bam_utils.write_position_fasta(input_bam_fp, temp_positions_fp, output_fasta_fp)

    os.remove(temp_positions_fp)


class BlastAnnotator(object):
    def __init__(self, annotations, database='GRCh38.d1.vd1.fa', max_target_seqs=5, max_hsps=5,
            rna_editing_identity_threshold=.95, rna_editing_coverage_threshold=.9):
        self.annotations = annotations
        self.database = database
        self.max_target_seqs = max_target_seqs
        self.max_hsps = max_hsps

        self.out_fields = get_required_out_fields(annotations)

        self.rna_editing_identity_threshold = rna_editing_identity_threshold
        self.rna_editing_coverage_threshold = rna_editing_coverage_threshold

    def blastn_fasta(self, input_fasta):
        """Blastn the given fasta and collect results for each sequence in input fasta"""
        blast_output = execute_blastn(input_fasta, self.database, out_fields=self.out_fields,
                max_target_seqs=self.max_target_seqs, max_hsps=self.max_hsps)
        output_dicts = parse_blast_output(blast_output, self.out_fields)

        sequence_to_results = defaultdict(list)
        for d in output_dicts:
            if 'qseqid' in d:
                sequence_to_results[d['qseqid']].append(d)

        return sequence_to_results

    def get_rna_editing_annotations(self, position_to_read_results):
        """Returns dict {position: %passing}"""
        position_to_percent_passing = {}
        for (chrom, pos), read_to_result_dicts in position_to_read_results.items():
            count = 0
            for result_dicts in read_to_result_dicts.values():
                if is_positive_rna_count(chrom, pos, result_dicts,
                        identity_threshold=self.rna_editing_identity_threshold,
                        coverage_threshold=self.rna_editing_coverage_threshold):
                    count += 1

            position_to_percent_passing[(chrom, pos)] = count / max(1, len(read_to_result_dicts))

        return position_to_percent_passing

    def get_rna_editing_blast_annotations(self, input_fasta, chrom_regex=r'^(.*):.*$',
            pos_regex=r'^.*:(.*)\|.*$', read_regex=r'^.*\|(.*)$'):
        """Collect by value in input fasta sequence identifier. Identifier is seperated by |

        i.e. if a sequence id is chr1:12345|read1, then the returned dictionary will
        look something like - {'chr1:12345': {read1: [{blastn parsed result}, ...], ...}, ...}
        """
        sequence_to_results = self.blastn_fasta(input_fasta)

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

    def get_blast_annotations_for_bam(self, input_bam_fp, position_tups):
        """Get annotations for the given positions based on reads in the given bam.

        input_bam_fp - filepath to bam with reads covering positions in the positions file
        position_tups - positions to recieve annotations. format - [(chrom, pos), ...]

        if position is not present in bam, then it will not be returned in annotations

        Returns: position_to_annotations, headers
            {(chrom, pos): [annotation1, annotation2, annotation3, ...]}, [header1, 
                    header2, header3, ...]
        """
        temp_fasta_fp = 'temp.query.fa'
        prepare_input_files(input_bam_fp, temp_fasta_fp, position_tups)

        annotations_dict = defaultdict(list)
        headers = []
        if 'rna_editing' in self.annotations:
            position_to_percent_passing = self.get_rna_editing_blast_annotations(temp_fasta_fp)
            # add positions that were missing for whatever reason
            for (chrom, pos) in position_tups:
                if (chrom, int(pos)) not in position_to_percent_passing:
                    position_to_percent_passing[(chrom, int(pos))] = '.'

            headers += ['BLAST_RNA_EDITING_%_PASSING']
            for (chrom, pos), value in position_to_percent_passing.items():
                annotations_dict[(chrom, str(pos))].append(value)

        os.remove(temp_fasta_fp)

        return annotations_dict, headers
