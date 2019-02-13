import os
import subprocess


GENE_TO_PRIMARY_TRANSCRIPT_FP = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        'data/gene_to_primary_transcript.tsv')
GENE_TO_PRIMARY_TRANSCRIPT = {l.strip().split('\t')[0]:l.strip().split('\t')[1] for l in open(GENE_TO_PRIMARY_TRANSCRIPT_FP)}

def get_gene_from_transvar_lines(transvar_lines):
    for line in transvar_lines:
        pieces = line.strip().split('\t')
        if len(pieces) > 6:
            gene = pieces[2]

            if gene in GENE_TO_PRIMARY_TRANSCRIPT:
                return gene

    return ''

def parse_for_ensembl_transcript(transvar_output, ensembl_transcript=None, use_primary=True):
    """Returns transcript, gene, strand, region, and info associated primary transcript or given
    transcript id.
    
    If transcript is not found, it returns first ensembl transcript in output.
    
    If no Ensembl transcript is found, returns first transcript it sees.
    
    Otherwise, returns (., ., .)"""


    lines = [l for l in transvar_output.split('\n')
            if 'input' not in l]

    main_gene = get_gene_from_transvar_lines(lines)
    ensembl_transcript = GENE_TO_PRIMARY_TRANSCRIPT.get(main_gene, '')

    first_transcript = None
    first_ensembl = None
    primary_transcript = None
    for i, line in enumerate(lines):
        pieces = line.strip().split('\t')
         # check for valid output
        if len(pieces) > 6:
            transcript = pieces[1].split(' ')[0]
            gene = pieces[2]
            strand = pieces[3]
            region = pieces[5]
            info = pieces[6]

            # check for first transcript
            if i == 0:
                first_transcript = (transcript, gene, strand, region, info)
            # check for ensembl
            if 'source=Ensembl' in info and first_ensembl is None:
                first_ensembl = (transcript, gene, strand, region, info)
            # check for primary transcript
            is_primary = transcript.lower().split('.')[0] == ensembl_transcript.lower().split('.')[0]
            if 'source=Ensembl' in info and is_primary:
                return transcript, gene, strand, region, info
            elif is_primary:
                primary_transcript = (transcript, gene, strand, region, info)

    if primary_transcript is not None:
        return primary_transcript
    if first_ensembl is not None:
        return first_ensembl
    if first_transcript is not None:
        return first_transcript
    return '.', '.', '.', '.', '.'

def get_transcript_gene_strand_region_info_tup(chrom, position, ensembl_transcript=None,
        use_primary=True):
    """Returns transcript, gene, strand, region, and info for the given position"""
    tool_args = ['transvar', 'ganno', '--ensembl', '--gencode', '--ucsc', '--refseq',
            '-i', f'{chrom}:g.{position}']
    result = subprocess.check_output(tool_args).decode('utf-8')

    if ensembl_transcript is None:
        return parse_for_ensembl_transcript(result, use_primary=use_primary)

    return parse_for_ensembl_transcript(result, ensembl_transcript=ensembl_transcript, use_primary=False)
