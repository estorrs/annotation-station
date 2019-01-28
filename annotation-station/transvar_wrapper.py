import subprocess

import subprocess

def parse_transvar_output(transvar_output):
    """Returns first transcript in transvar output"""
    for line in transvar_output.split('\n'):
        if 'input' not in line:
            pieces = line.strip().split('\t')

            if len(pieces) < 6:
                return '.', '.', '.'

            gene = pieces[2]
            strand = pieces[3]
            region = pieces[5]
            return gene, strand, region

    return '.', '.', '.'



def get_gene_strand_region_tup(chrom, position):
    """Returns gene, strand, and region for the given position"""
    tool_args = ['transvar', 'ganno', '--ccds',
            '-i', f'{chrom}:g.{position}']
    result = subprocess.check_output(tool_args).decode('utf-8')

    return parse_transvar_output(result)
