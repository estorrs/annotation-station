import argparse
import os
import subprocess

import transvar_wrapper

parser = argparse.ArgumentParser()

parser.add_argument('input_file', type=str,
        help='input file')

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

def check_transvar_setup():
    """Will set up transvar if needed"""
    try:
        result = transvar_wrapper.get_gene_strand_region_tup('chr1', '12345')
    except subprocess.CalledProcessError:
        print('Setting up transvar')
        tool_args = ['transvar', 'config', '--download_ref', '--refversion', 'hg38']
        subprocess.check_output(tool_args)
    
        tool_args = ['transvar', 'config', '--download_anno', '--refversion', 'hg38']
        subprocess.check_output(tool_args)
        print('finished transvar setup')


def annotate_tsv(input_fp, output_fp, input_header=False):
    f = open(input_fp)
    out_f = open(output_fp, 'w')
    
    if input_header:
        out_f.write(f.readline()[:-1] + '\tGENE\tSTRAND\tPOSITION\n')
    
    for line in f:
        pieces = line.strip().split('\t', 2)
        chrom = pieces[0]
        pos = pieces[1]

        gene, strand, region = transvar_wrapper.get_gene_strand_region_tup(chrom, pos)
        
        out_f.write(line[:-1] + f'\t{gene}\t{strand}\t{region}\n')

    f.close()
    out_f.close()

def main():
    check_arguments()

    check_transvar_setup()

    annotate_tsv(args.input_file, args.output, input_header=args.input_header)

if __name__ == '__main__':
    main()
