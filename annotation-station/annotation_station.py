import argparse
import os
import subprocess

from transvar_wrapper import TransvarAnnotator
from repeats import RepeatAnnotator

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
annotation_group.add_argument('--annotate-blast', action='store_true',
        help='If present, annotations for BLAST will be done. \
Added fields will include HI')

parser.add_argument('--primary-transcripts', type=str,
        help='A .tsv file with genes in first column and ensembl transcript in second column. \
Only necessary if --annotate-transvar flag is present')
parser.add_argument('--repeats-table', type=str,
        help='A .tsv file generated with ucsc table browser - repeats. \
Only necessary if --annotate-repeats flag is present')
parser.add_argument('--input-header', action='store_true',
        help='Whether input tsv file has header or not')
parser.add_argument('--output', type=str,
        default='output.tsv', help='output fp')
parser.add_argument('--input-type', type=str,
        help='Type of input file. Options are tsv and json.')
parser.add_argument('--reference-version', type=str,
        default='hg38', help='Reference version to use for annotations')

args = parser.parse_args()

DEFAULT_REPEATS_TABLE = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'data/repeats_table.grch38.tsv')

def check_arguments():
    if args.input_type is None:
        raise ValueError('Must specify an input type')

    if not args.annotate_transvar and not args.annotate_repeats:
        raise ValueError('Must specify a type of annotation to perform.')

    if args.annotate_transvar and args.primary_transcripts is None:
        raise ValueError('Must specify a primary transcripts file with --primary-transcripts flag \
if using --annotate-transvar.')

#     if args.annotate_repeats and args.repeats_table is None:
#         raise ValueError('Must specify a repeats file with --repeats-table flag \
# if using --annotate-repeats. File can be downloaded with ucsc table browser (repeats).')

def check_transvar_setup(transvar_annotator, reference_version='hg38'):
    """Will set up transvar if needed"""
    try:
        result = transvar_annotator.get_transcript_gene_strand_region_info_tup('chr1', '12345',
                reference_version=reference_version)
    except subprocess.CalledProcessError:
        print('Setting up transvar')
        tool_args = ['transvar', 'config', '--download_ref', '--refversion', reference_version]
        subprocess.check_output(tool_args)
    
        tool_args = ['transvar', 'config', '--download_anno', '--refversion', reference_version]
        subprocess.check_output(tool_args)
        print('finished transvar setup')

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

def annotate_transvar_tsv(transvar_annotator, fp, input_header=False, reference_version='hg38'):
    """Annotate transvar tsv"""
    out_lines = []

    f = open(fp)

    if input_header:
        out_lines.append(f.readline()[:-1] + '\tPRIMARY_TRANSCRIPT\tGENE\tSTRAND\tREGION\tNON_VERBOSE_REGION\tINFO')

    for line in f:
        pieces = line.strip().split('\t', 2)
        chrom = pieces [0]
        pos = pieces[1]

        transcript, gene, strand, region, info = transvar_annotator.get_transcript_gene_strand_region_info_tup(
                chrom, pos, reference_version=reference_version)
        simplified_region = get_simplified_region(region)

        out_lines.append(line[:-1] + f'\t{transcript}\t{gene}\t{strand}\t{region}\t{simplified_region}\t{info}')
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

def main():
    check_arguments()

    # create our output file
    f = open(args.input_file)
    out_f = open(args.output, 'w')
    out_f.write(f.read())
    f.close()
    out_f.close()

    if args.annotate_transvar:
        ta = TransvarAnnotator(args.primary_transcripts)
        check_transvar_setup(ta, reference_version=args.reference_version)
        annotate_transvar_tsv(ta, args.output, input_header=args.input_header,
                reference_version=args.reference_version)
    if args.annotate_repeats:
        if args.repeats_table is None:
            ra = RepeatAnnotator(DEFAULT_REPEATS_TABLE)
        else:
            ra = RepeatAnnotator(args.repeats_table)
        annotate_repeats_tsv(ra, args.output, input_header=args.input_header)
#     if args.annotate_blast:
#         annotate_blast_tsv(args.output)

if __name__ == '__main__':
    main()
