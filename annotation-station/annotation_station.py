import argparse
import os
import subprocess

#from blast import BlastAnnotator
from blat import BlatAnnotator
from repeats import RepeatAnnotator
from transvar_wrapper import TransvarAnnotator

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
# annotation_group.add_argument('--annotate-blast', action='store_true',
#         help='If present, annotations for BLAST will be done. \
# Added fields will include BLAST_RNA_EDITING_%_PASSING')
annotation_group.add_argument('--annotate-blat', action='store_true',
        help='If present, annotations for BLAT will be done. \
Added fields will include BLAT_RNA_EDITING_%_PASSING')

# transvar specific
parser.add_argument('--primary-transcripts', type=str,
        help='A .tsv file with genes in first column and ensembl transcript in second column. \
Only used if --annotate-transvar flag is present')

# repeats specific
parser.add_argument('--repeats-table', type=str,
        help='A .tsv file generated with ucsc table browser - repeats. \
Only used if --annotate-repeats flag is present')

# blat specific
parser.add_argument('--blat-database', type=str,
        help='Database to use with blat.')
parser.add_argument('--blat-input-bam', type=str,
        help='Input bam containing reads to use for blat annotation. \
Required if --annotate-blat is used.')
parser.add_argument('--rna-editing-percent-threshold', type=float,
        default=.95, help='Percent identity threshold to use when calling a positive blat rna \
editing read.')

# # blast specific
# parser.add_argument('--blast-database', type=str,
#         help='Database to use with blast. Database must be installed on host')
# parser.add_argument('--blast-input-bam', type=str,
#         help='Input bam containing reads to use for blast annotation. \
# Required if --annotate-blast is used.')
# parser.add_argument('--rna-editing-identity-threshold', type=float,
#         default=.95, help='Percent identity threshold to use when calling a positive blast rna \
# editing read.')
# parser.add_argument('--rna-editing-coverage-threshold', type=float,
#         default=.90, help='Percent sequence coverage threshold to use when calling a positive \
# blast rna editing read.')

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
# DEFAULT_GRCH38_BLAST_DATABASE = 'GRCh38.d1.vd1.fa'
# DEFAULT_GRCH37_BLAST_DATABASE = 'ucsc.hg19.fa'

def check_arguments():
    if args.input_type is None:
        raise ValueError('Must specify an input type')

#     if args.annotate_blast and args.blast_input_bam is None:
#         raise ValueError('Must specify input bam to use for blast annotations')

#     if not args.annotate_transvar and not args.annotate_repeats:
#         raise ValueError('Must specify a type of annotation to perform.')
# 
#     if args.annotate_transvar and args.primary_transcripts is None:
#         raise ValueError('Must specify a primary transcripts file with --primary-transcripts flag \
# if using --annotate-transvar.')

#     if args.annotate_repeats and args.repeats_table is None:
#         raise ValueError('Must specify a repeats file with --repeats-table flag \
# if using --annotate-repeats. File can be downloaded with ucsc table browser (repeats).')

def get_default_repeat_table(reference_version):
    """Returns default repeat table fp for given reference"""
    if reference_version == 'hg38' or reference_version == 'grch38':
        return DEFAULT_GRCH38_REPEATS_TABLE
    if reference_version == 'hg19' or reference_version == 'hg37' or reference_version == 'grch37':
        return DEFAULT_GRCH37_REPEATS_TABLE
    raise ValueError('Incompatible reference version for built in repeats table')

# def get_default_blast_database(reference_version):
#     """Return default blast database for given reference"""
#     if reference_version == 'hg38' or reference_version == 'grch38':
#         return DEFAULT_GRCH38_BLAST_DATABASE
#     if reference_version == 'hg19' or reference_version == 'hg37' or reference_version == 'grch37':
#         return DEFAULT_GRCH37_BLAST_DATABASE
#     raise ValueError('Incompatible reference version for built in blast database')

def transvar_setup(transvar_annotator, reference_version='hg38', reference_fasta=None):
    """Will set up transvar if needed"""
    try:
        result = transvar_annotator.get_transcript_gene_strand_region_info_tup('chr1', '12345',
                reference_version=reference_version)
    except subprocess.CalledProcessError:
        print('Setting up transvar')
        if reference_fasta is None:
            tool_args = ['transvar', 'config', '--download_ref', '--refversion', reference_version]
        else:
            tool_args = ['transvar', 'config', '-k', 'reference',
                    '-v', reference_fasta,
                    '--refversion', reference_version]
        subprocess.check_output(tool_args)
    
        tool_args = ['transvar', 'config', '--download_anno', '--refversion', reference_version]
        subprocess.check_output(tool_args)
        print('finished transvar setup')

# def check_transvar_setup(transvar_annotator, reference_version='hg38'):
#     """Will set up transvar if needed"""
#     try:
#         result = transvar_annotator.get_transcript_gene_strand_region_info_tup('chr1', '12345',
#                 reference_version=reference_version)
#     except subprocess.CalledProcessError:
#         print('Setting up transvar')
#         tool_args = ['transvar', 'config', '--download_ref', '--refversion', reference_version]
#         subprocess.check_output(tool_args)
#     
#         tool_args = ['transvar', 'config', '--download_anno', '--refversion', reference_version]
#         subprocess.check_output(tool_args)
#         print('finished transvar setup')

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

def annotate_blat_tsv(blat_annotator, fp, input_bam, input_header=False):
    out_lines = []
    f = open(fp)
    if input_header:
        f.readline()
    chrom_pos_tups = []
    for line in f:
        chrom, pos = line.strip().split('\t', 2)[:2]
        chrom_pos_tups.append((chrom, pos))
    f.close()

    blat_annotations_dict, headers = blat_annotator.get_blat_annotations_for_bam(
            input_bam, chrom_pos_tups)

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

    if args.annotate_transvar:
        if args.primary_transcripts is None:
            ta = TransvarAnnotator(DEFAULT_GENE_TO_PRIMARY_TRANSCRIPT)
        else:
            ta = TransvarAnnotator(args.primary_transcripts)
        check_transvar_setup(ta, reference_version=args.reference_version,
                reference_fasta=args.reference_fasta)
        annotate_transvar_tsv(ta, args.output, input_header=args.input_header,
                reference_version=args.reference_version)

    if args.annotate_repeats:
        if args.repeats_table is None:
            ra = RepeatAnnotator(get_default_repeat_table(args.reference_version))
        else:
            ra = RepeatAnnotator(args.repeats_table)
        annotate_repeats_tsv(ra, args.output, input_header=args.input_header)

    if args.annotate_blat:
        ba = BlatAnnotator(['rna_editing'],
                database=args.blat_database,
                rna_editing_percent_threshold=args.rna_editing_percent_threshold)
        annotate_blat_tsv(ba, args.output, args.blast_input_bam,
                input_header=args.input_header)

if __name__ == '__main__':
    main()
