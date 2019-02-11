import argparse
import os
import subprocess

import transvar_wrapper

parser = argparse.ArgumentParser()

parser.add_argument('input_file', type=str,
        help='input file')

annotation_group = parser.add_argument_group('argument_group') 
annotation_group.add_argument('--annotate-transvar', action='store_true',
        help='If present, transvar annotations will be added to input file. \
If --ensembl-transcript is also provided, the annotation for that transcript will be done. \
Added fields will be TRANSCRIPT, GENE, STRAND, REGION, and INFO.')
annotation_group.add_argument('--annotate-alu', action='store_true',
        help='If present, ALU annotation will be done. \
Boolean IS_ALU and ALU_SUBTYPE fields will be added.')
annotation_group.add_argument('--annotate-blast', action='store_true',
        help='If present, annotations for BLAST will be done. \
Added fields will include HI')

parser.add_argument('--input-header', action='store_true',
        help='Whether input tsv file has header or not')
parser.add_argument('--output', type=str,
        default='output.tsv', help='output fp')
parser.add_argument('--input-type', type=str,
        help='Type of input file. Options are tsv and json.')
parser.add_argument('--ensembl-transcript', type=str,
        default=None, help='Ensembl transcript to use with transvar annotations')

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

def annotate_transvar_tsv(fp, input_header=False, ensembl_transcript=None):
    """Annotate transvar tsv"""
    out_lines = []

    f = open(fp)

    if input_header:
        out_lines.append(f.readline()[:-1] + '\tTRANSCRIPT\tGENE\tSTRAND\tREGION\tINFO\n')

    for line in f:
        pieces = line.strip().split('\t', 2)
        chrom = pieces [0]
        pos = pieces[1]

        transcript, gene, strand, region, info = transvar_wrapper.get_transcript_gene_strand_region_info_tup(
                chrom, pos, ensembl_transcript=ensembl_transcript)

        out_lines.append(line[:-1] + f'\t{transcript}\t{gene}\t{strand}\t{region}\t{info}\n')
    f.close()

    output_str = '\n'.join(out_lines) + '\n'
    # write over old file
    f = open(fp, 'w')
    f.write(output_str)
    f.close()

# def annotate_tsv(input_fp, output_fp, input_header=False):
#     f = open(input_fp)
#     out_f = open(output_fp, 'w')
#     
#     if input_header:
#         out_f.write(f.readline()[:-1] + '\tGENE\tSTRAND\tPOSITION\n')
#     
#     for line in f:
#         pieces = line.strip().split('\t', 2)
#         chrom = pieces[0]
#         pos = pieces[1]
# 
#         gene, strand, region = transvar_wrapper.get_gene_strand_region_tup(chrom, pos)
#         
#         out_f.write(line[:-1] + f'\t{gene}\t{strand}\t{region}\n')
# 
#     f.close()
#     out_f.close()

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
        annotate_transvar_tsv(args.output, input_header=args.input_header,
            ensembl_transcript=args.ensembl_transcript)
#     if args.annotate_alu:
#         annotate_alu_tsv(args.output)
#     if args.annotate_blast:
#         annotate_blast_tsv(args.output)

if __name__ == '__main__':
    main()
