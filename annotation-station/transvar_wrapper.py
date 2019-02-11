import subprocess


def parse_for_ensembl_transcript(transvar_output, ensembl_transcript):
    """Returns transcript, gene, strand, region, and info associated with the given ensembl
    transcript id.
    
    If transcript is not found, it returns first ensembl transcript in output.
    
    If no Ensembl transcript is found, returns first transcript it sees.
    
    Otherwise, returns (., ., .)"""
    lines = [l for l in transvar_output.split('\n')
            if 'input' not in l]

    first_transcript = None
    first_ensembl = None
    primary_transcript = None
    for i, line in enumerate(lines):
        pieces = line.strip().split('\t')
         # check for valid output
        if len(pieces) > 6:
            transcript = pieces[1]
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
            if 'source=Ensembl' in info and transcript.lower() == ensembl_transcript.lower().split('.')[0]:
                return transcript, gene, strand, region, info
            elif transcript.lower() == ensembl_transcript.lower().split('.')[0]:
                primary_transcript = (transcript, gene, strand, region, info)

    if primary_transcript is not None:
        return primary_transcript
    if first_ensemble is not None:
        return first_ensembl
    if first_transcript is not None:
        return first_transcript
    return '.', '.', '.', '.', '.'

def get_transcript_gene_strand_region_info_tup(chrom, position, ensembl_transcript=None):
    """Returns transcript, gene, strand, region, and info for the given position"""
    tool_args = ['transvar', 'ganno', '--ensembl', '--gencode', '--ucsc', '--refseq',
            '-i', f'{chrom}:g.{position}']
    result = subprocess.check_output(tool_args).decode('utf-8')

    if ensembl_transcript is None:
        return parse_for_ensembl_transcript(result, '')
    return parse_for_ensembl_transcript(result, ensembl_transcript)
