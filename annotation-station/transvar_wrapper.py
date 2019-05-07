import os
import subprocess

class TransvarAnnotator(object):
    def __init__(self, gene_to_primary_transcript_fp):
        self.gene_to_primary_transcript = {l.strip().split('\t')[0]:l.strip().split('\t')[1]
                for l in open(gene_to_primary_transcript_fp)}

    def get_gene_from_transvar_lines(self, transvar_lines):
        for line in transvar_lines:
            pieces = line.strip().split('\t')
            if len(pieces) > 6:
                gene = pieces[2]

                if gene in self.gene_to_primary_transcript:
                    return gene

        return ''

    def parse_for_ensembl_transcript(self, transvar_output, ensembl_transcript=None, use_primary=True):
        """Returns transcript, gene, strand, region, and info associated primary transcript or given
        transcript id.
    
        If transcript is not found, it returns first ensembl transcript in output.
    
        If no Ensembl transcript is found, returns first transcript it sees.
    
        Otherwise, returns (., ., .)"""


        lines = [l for l in transvar_output.split('\n')
                if 'input' not in l]

        main_gene = self.get_gene_from_transvar_lines(lines)
        ensembl_transcript = self.gene_to_primary_transcript.get(main_gene, '')

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
                coordinates = pieces[4]
                region = pieces[5]
                info = pieces[6]

                # check for first transcript
                if i == 0:
                    first_transcript = (transcript, gene, strand, coordinates, region, info)
                # check for ensembl
                if 'source=Ensembl' in info and first_ensembl is None:
                    first_ensembl = (transcript, gene, strand, coordinates, region, info)
                # check for primary transcript
                is_primary = transcript.lower().split('.')[0] == ensembl_transcript.lower().split('.')[0]
                if 'source=Ensembl' in info and is_primary:
                    return transcript, gene, strand, coordinates, region, info
                elif is_primary:
                    primary_transcript = (transcript, gene, strand, coordinates, region, info)

        if primary_transcript is not None:
            return primary_transcript
        if first_ensembl is not None:
            return first_ensembl
        if first_transcript is not None:
            return first_transcript
        return '.', '.', '.', '.', '.', '.'

    def get_transcript_gene_strand_region_info_tup(self, chrom, position, ensembl_transcript=None,
            use_primary=True, reference_version='hg38', ref_base=None, alt_base=None):
        """Returns transcript, gene, strand, region, and info for the given position"""
        if ref_base is not None and alt_base is not None:
            tool_args = ['transvar', 'ganno', '--ensembl', '--gencode', '--ucsc', '--refseq',
                    '-i', f'{chrom}:g.{position}{ref_base.upper()}>{alt_base.upper()}',
                    '--refversion', reference_version]
        else:
            tool_args = ['transvar', 'ganno', '--ensembl', '--gencode', '--ucsc', '--refseq',
                    '-i', f'{chrom}:g.{position}',
                    '--refversion', reference_version]
        result = subprocess.check_output(tool_args).decode('utf-8')

        if ensembl_transcript is None:
            return self.parse_for_ensembl_transcript(result, use_primary=use_primary)

        return self.parse_for_ensembl_transcript(result, ensembl_transcript=ensembl_transcript, use_primary=False)
