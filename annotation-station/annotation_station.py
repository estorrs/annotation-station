import argparse
import os
import subprocess

import transvar_wrapper
import repeats

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

parser.add_argument('--input-header', action='store_true',
        help='Whether input tsv file has header or not')
parser.add_argument('--output', type=str,
        default='output.tsv', help='output fp')
parser.add_argument('--input-type', type=str,
        help='Type of input file. Options are tsv and json.')

args = parser.parse_args()

def check_arguments():
    if args.input_type is None:
        raise ValueError('Must specify an input type')
    if args.annotate_transvar is None and args.annotate_alu is None and args.annotate_blast is None:
        raise ValueError('Must specify a type of annotation to perform.')

def check_transvar_setup():
    """Will set up transvar if needed"""
    try:
        result = transvar_wrapper.get_transcript_gene_strand_region_info_tup('chr1', '12345')
    except subprocess.CalledProcessError:
        print('Setting up transvar')
        tool_args = ['transvar', 'config', '--download_ref', '--refversion', 'hg38']
        subprocess.check_output(tool_args)
    
        tool_args = ['transvar', 'config', '--download_anno', '--refversion', 'hg38']
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

def annotate_transvar_tsv(fp, input_header=False):
    """Annotate transvar tsv"""
    out_lines = []

    f = open(fp)

    if input_header:
        out_lines.append(f.readline()[:-1] + '\tPRIMARY_TRANSCRIPT\tGENE\tSTRAND\tREGION\tNON_VERBOSE_REGION\tINFO')

    for line in f:
        pieces = line.strip().split('\t', 2)
        chrom = pieces [0]
        pos = pieces[1]

        transcript, gene, strand, region, info = transvar_wrapper.get_transcript_gene_strand_region_info_tup(
                chrom, pos)
        simplified_region = get_simplified_region(region)

        out_lines.append(line[:-1] + f'\t{transcript}\t{gene}\t{strand}\t{region}\t{simplified_region}\t{info}')
    f.close()

    output_str = '\n'.join(out_lines) + '\n'
    # write over old file
    f = open(fp, 'w')
    f.write(output_str)
    f.close()

def annotate_repeats_tsv(fp, input_header=False):
    """Annotate repeats tsv"""
    out_lines = []

    f = open(fp)

    if input_header:
        out_lines.append(f.readline()[:-1] + '\tREPEAT_NAME\tREPEAT_CLASS\tREPEAT_FAMILY')

    for line in f:
        pieces = line.strip().split('\t', 2)
        chrom = pieces [0]
        pos = pieces[1]

        repeat = repeats.get_repeat_by_position(chrom, pos)

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
        check_transvar_setup()
        annotate_transvar_tsv(args.output, input_header=args.input_header)
    if args.annotate_repeats:
        annotate_repeats_tsv(args.output, input_header=args.input_header)
#     if args.annotate_blast:
#         annotate_blast_tsv(args.output)

if __name__ == '__main__':
    main()
