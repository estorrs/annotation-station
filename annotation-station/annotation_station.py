import argparse
import logging
import os
import subprocess

import bam_utils
#from blast import BlastAnnotator
from blat import BlatAnnotator
from repeats import RepeatAnnotator
from transvar_wrapper import TransvarAnnotator

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('input_file', type=str,
        help='input file')

annotation_group = parser.add_argument_group('argument_group') 
annotation_group.add_argument('--annotate-transvar', action='store_true',
        help='If present, transvar annotations will be added to input file. \
Added fields will be PRIMARY_TRANSCRIPT, GENE, STRAND, REGION, NON_VERBOSE_REGION, and INFO.')
annotation_group.add_argument('--annotate-repeats', action='store_true',
        help='If present, annotation for repeats will be done. i.e. for finding ALUs and stuff. \
Added fields will be REPEAT_NAME, REPEAT_CLASS, and REPEAT_FAMILY. \
Entry will be . if position is not a repeat.')
annotation_group.add_argument('--annotate-blat', action='store_true',
        help='If present, annotations for BLAT will be done. \
Added fields will include BLAT_RNA_EDITING_%_PASSING')

# transvar specific
parser.add_argument('--primary-transcripts', type=str,
        help='A .tsv file with genes in first column and ensembl transcript in second column. \
Only used if --annotate-transvar flag is present')
parser.add_argument('--with-base-change', action='store_true',
        help='Use base ref/alt base change with annotations. The third column in input file must \
be the reference base, and the fourth column must be the alternative base')

# repeats specific
parser.add_argument('--repeats-table', type=str,
        help='A .tsv file generated with ucsc table browser - repeats. \
Only used if --annotate-repeats flag is present')

# blat specific
parser.add_argument('--blat-input-bam', type=str,
        help='Input bam containing reads to use for blat annotation. \
Required if --annotate-blat is used.')
parser.add_argument('--rna-editing-percent-threshold', type=float,
        default=.95, help='Percent identity threshold to use when calling a positive blat rna \
editing read.')


parser.add_argument('--input-header', action='store_true',
        help='Whether input tsv file has header or not')
parser.add_argument('--input-type', type=str,
        help='Type of input file. Options are tsv and json.')
parser.add_argument('--reference-version', type=str,
        default='hg38', help='Reference version to use for annotations. \
Important for repeats and transvar')
parser.add_argument('--reference-fasta', type=str,
        help='Reference fasta to use for annotations with blat and transvar')
parser.add_argument('--output', type=str,
        default='output.tsv', help='output fp')

args = parser.parse_args()


# defaults
DEFAULT_GRCH38_REPEATS_TABLE = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'data/repeats/repeats_table.grch38.tsv')
DEFAULT_GRCH37_REPEATS_TABLE = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'data/repeats/repeats_table.hg19.tsv')
DEFAULT_GENE_TO_PRIMARY_TRANSCRIPT = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'data/transcripts/gene_to_primary_transcript.tsv')
DEFAULT_GENE_TO_PRIMARY_TRANSCRIPT = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'data/transcripts/gene_to_primary_transcript.tsv')

def check_arguments():
    if args.input_type is None:
        raise ValueError('Must specify an input type')

def get_default_repeat_table(reference_version):
    """Returns default repeat table fp for given reference"""
    if reference_version.lower() == 'hg38' or reference_version.lower() == 'grch38':
        return DEFAULT_GRCH38_REPEATS_TABLE
    if reference_version.lower() == 'hg19' or reference_version.lower() == 'hg37' or reference_version.lower() == 'grch37':
        return DEFAULT_GRCH37_REPEATS_TABLE
    raise ValueError('Incompatible reference version for built in repeats table')

def check_transvar_setup(transvar_annotator, reference_version='hg38', reference_fasta=None):
    """Will set up transvar if needed"""
    try:
        result = transvar_annotator.get_transcript_gene_strand_region_info_tup('chr1', '12345',
                reference_version=reference_version)
    except subprocess.CalledProcessError:
        logging.info('Setting up transvar')
        if reference_fasta is None:
            tool_args = ['transvar', 'config', '--download_ref', '--refversion', reference_version]
        else:
            tool_args = ['transvar', 'config', '-k', 'reference',
                    '-v', reference_fasta,
                    '--refversion', reference_version]
        subprocess.check_output(tool_args)
    
        tool_args = ['transvar', 'config', '--download_anno', '--refversion', reference_version]
        subprocess.check_output(tool_args)
        logging.info('finished transvar setup')

def get_simplified_region(transvar_region):
    """Converts transvar region to a simplified region

    Possible values: CODING_EXONIC, NON_CODING_EXONIC, INTRONIC, 5UTR, 3UTR, INTERGENIC, .
    """
    if 'cds_in_exon' in transvar_region:
        return 'CODING_EXONIC'
    elif '5-UTR' in transvar_region:
        return '5UTR'
    elif '3-UTR' in transvar_region:
        return '3UTR'
    elif 'noncoding_exon' in transvar_region:
        return 'NON_CODING_EXONIC'
    elif 'intron' in transvar_region:
        return 'INTRONIC'
    elif 'intergenic' in transvar_region:
        return 'INTERGENIC'
    return '.'

def annotate_transvar_tsv(transvar_annotator, fp, input_header=False, reference_version='hg38',
            with_base_change=False):
    """Annotate transvar tsv"""
    out_lines = []

    f = open(fp)

    if input_header:
        out_lines.append(f.readline()[:-1] + '\tPRIMARY_TRANSCRIPT\tGENE\tSTRAND\tCOORDINATES\tREGION\tNON_VERBOSE_REGION\tINFO')

    for line in f:
        pieces = line.strip().split('\t')
        chrom = pieces[0]
        pos = pieces[1]

        if with_base_change:
            ref_base, alt_base = pieces[2], pieces[3]
            transcript, gene, strand, coordinates, region, info = transvar_annotator.get_transcript_gene_strand_region_info_tup(
                    chrom, pos, reference_version=reference_version, ref_base=ref_base,
                    alt_base=alt_base)
        else:
            transcript, gene, strand, coordinates, region, info = transvar_annotator.get_transcript_gene_strand_region_info_tup(
                    chrom, pos, reference_version=reference_version)

        simplified_region = get_simplified_region(region)

        out_lines.append(line[:-1] + f'\t{transcript}\t{gene}\t{strand}\t{coordinates}\t{region}\t{simplified_region}\t{info}')
    f.close()

    output_str = '\n'.join(out_lines) + '\n'
    # write over old file
    f = open(fp, 'w')
    f.write(output_str)
    f.close()

def annotate_repeats_tsv(repeat_annotator, fp, input_header=False):
    """Annotate repeats tsv"""
    out_lines = []

    f = open(fp)

    if input_header:
        out_lines.append(f.readline()[:-1] + '\tREPEAT_NAME\tREPEAT_CLASS\tREPEAT_FAMILY')

    for line in f:
        pieces = line.strip().split('\t', 2)
        chrom = pieces [0]
        pos = pieces[1]

        repeat = repeat_annotator.get_repeat_by_position(chrom, pos)

        if repeat is not None:
            repeat_name, repeat_class, repeat_family = repeat
        else:
            repeat_name, repeat_class, repeat_family = '.', '.', '.'

        out_lines.append(line[:-1] + f'\t{repeat_name}\t{repeat_class}\t{repeat_family}')
    f.close()

    output_str = '\n'.join(out_lines) + '\n'
    # write over old file
    f = open(fp, 'w')
    f.write(output_str)
    f.close()

def annotate_blat_tsv(blat_annotator, fp, input_bam, input_header=False):
    out_lines = []
    f = open(fp)
    if input_header:
        f.readline()
    chrom_pos_tups = []
    reference_bases = []
    for line in f:
        chrom, pos, base = line.strip().split('\t', 3)[:3]
        chrom_pos_tups.append((chrom, pos))
        reference_bases.append(base)
    f.close()

    # chunk into 5000 position pieces
    chunk_size = 1000
    chunked_chrom_pos_tups = []
    chunked_reference_bases = []
    prev = 0
    logging.info(f'starting processing of {len(chrom_pos_tups)} total positions')
    logging.info(f'using chunk size of {chunk_size}')
    for i in range(chunk_size, len(chrom_pos_tups) + chunk_size, chunk_size):
        chunked_chrom_pos_tups.append(chrom_pos_tups[prev:i])
        chunked_reference_bases.append(reference_bases[prev:i])
        prev = i
    
    blat_annotations_dict = {}
    headers = []
    for i, (chrom_pos_chunk, reference_bases_chunk) in enumerate(zip(
            chunked_chrom_pos_tups, chunked_reference_bases)):
        if chrom_pos_chunk and reference_bases_chunk:
            logging.info(f'processing chunk {i + 1} of {len(chunked_chrom_pos_tups)}')
            d, h = blat_annotator.get_blat_annotations_for_bam(input_bam, chrom_pos_chunk,
                    reference_bases=reference_bases_chunk)
            blat_annotations_dict.update(d)
            headers = h

    f = open(fp)
    if input_header:
        out_lines.append(f.readline()[:-1] + '\t' + '\t'.join(headers))

    for i, line in enumerate(f):
        chrom, pos = chrom_pos_tups[i]
        annotations = blat_annotations_dict[(chrom, pos)]
        out_lines.append(line[:-1] + '\t' + '\t'.join([str(a) for a in annotations]))
    f.close()

    output_str = '\n'.join(out_lines) + '\n'
    # write over old file
    f = open(fp, 'w')
    f.write(output_str)
    f.close()

def main():
    check_arguments()

    # create our output file
    f = open(args.input_file)
    out_f = open(args.output, 'w')
    out_f.write(f.read())
    f.close()
    out_f.close()

    # index reference if it's there
    if args.reference_fasta is not None:
        bam_utils.index_reference(args.reference_fasta)

    if args.annotate_transvar:
        logging.info('Beginning transvar annotations')
        if args.primary_transcripts is None:
            ta = TransvarAnnotator(DEFAULT_GENE_TO_PRIMARY_TRANSCRIPT)
        else:
            ta = TransvarAnnotator(args.primary_transcripts)
        check_transvar_setup(ta, reference_version=args.reference_version,
                reference_fasta=args.reference_fasta)
        annotate_transvar_tsv(ta, args.output, input_header=args.input_header,
                reference_version=args.reference_version, with_base_change=args.with_base_change)

    if args.annotate_repeats:
        logging.info('Begginging repeat annotations')
        if args.repeats_table is None:
            ra = RepeatAnnotator(get_default_repeat_table(args.reference_version))
        else:
            ra = RepeatAnnotator(args.repeats_table)
        annotate_repeats_tsv(ra, args.output, input_header=args.input_header)

    if args.annotate_blat:
        logging.info('Beginning blat annotations')
        ba = BlatAnnotator(['rna_editing'],
                database=args.reference_fasta,
                rna_editing_percent_threshold=args.rna_editing_percent_threshold)
        annotate_blat_tsv(ba, args.output, args.blat_input_bam,
                input_header=args.input_header)

if __name__ == '__main__':
    main()
