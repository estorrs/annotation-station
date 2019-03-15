import os
import subprocess

def index_reference(reference_fasta_fp):
    """index the given reference if it doesnt exist"""
    if not os.path.isfile(reference_fasta_fp + '.fai'):
        tool_args = ['samtools', 'faidx', reference_fasta_fp]
        print(subprocess.check_output(tool_args).decode('utf-8'))

def index_bam(bam_fp):
    """index the given bam."""
    tool_args = ['samtools', 'index', bam_fp]
    print(subprocess.check_output(tool_args).decode('utf-8'))

def filter_bam_by_positions(bam_fp, positions_fp, output_fp, threads=1):
    """run bam filter step"""
    tool_args = ['samtools', 'view', '-h',
            '-@', str(threads),
            '-L', positions_fp,
            '-o', output_fp,
            bam_fp]

    print(subprocess.check_output(tool_args).decode('utf-8'))

class ReadCollection(object):
    def __init__(self, position_tups):
        """
        position tups - [(chrom, pos), ...]
        """
        self.positions = sorted(position_tups, key=lambda x: (x[0], x[1]))
        
        self.chrom_to_positions = {c:[] for c, _ in position_tups}
        for chrom, pos in position_tups:
            self.chrom_to_positions[chrom].append((chrom, pos))
        for chrom, ps in self.chrom_to_positions.items():
            self.chrom_to_positions[chrom] = sorted(ps, key=lambda x: (x[0], x[1]))

        self.position_to_reads = {p:[] for p in position_tups}

    def put_read(self, chrom, start, sequence):
        end = start + len(sequence)

        for chrom, pos in self.chrom_to_positions[chrom]:
            if pos > start and pos < end:
                base = sequence[pos - start]
                self.position_to_reads[(chrom, pos)].append((chrom, start, sequence, base))

            if pos > end:
                break

    def get_reads(self, chrom, pos):
        return self.position_to_reads.get((chrom, pos), [])


def get_reads_to_sequences_from_fasta(input_fasta_fp):
    f = open(input_fasta_fp)

    active_seq_id = None
    active_seq = ''
    reads_to_sequences = {}
    for line in f:
        line = line.strip()
        if line[0] == '>':
            if active_seq_id is not None:
                reads_to_sequences[active_seq_id] = active_seq
            active_seq_id = line[1:]
        else:
            active_seq += line
    f.close()

    return reads_to_sequences

def write_position_fasta(input_bam_fp, positions_fp, output_fasta_fp, max_depth=1000):
    """Writes a fasta with the given positions and bam.

    Will also return a dict mapping reads to their sequence"""
    tool_args = ['samtools', 'view',
            '-L', positions_fp,
             input_bam_fp]
    ps_1 = subprocess.Popen(tool_args, stdout=subprocess.PIPE)
    output = subprocess.check_output(('cut', '-f', '3,4,10'), stdin=ps_1.stdout).decode('utf-8')
    ps_1.wait()
    
    read_tups = []
    for line in output.split('\n'):
        # line is chrom\tstart\tsequence\n
        pieces = line.split('\t')
        if len(pieces) == 3:
            read_tups.append((pieces[0], int(pieces[1]), pieces[2]))

    # grab positions from file
    f = open(positions_fp)
    positions = []
    for line in f:
        # line is chrom\tpos....\n
        pieces = line.strip().split('\t')
        positions.append((pieces[0], int(pieces[1])))
    f.close()

    # create read collection
    rc = ReadCollection(positions)
    for chrom, start, seq in read_tups:
        rc.put_read(chrom, start, seq)
    
    f = open(output_fasta_fp, 'w')
    reads_to_data = {}
    for chrom, pos in positions:
        reads = rc.get_reads(chrom, pos)
        for i, (_, _, seq, base) in enumerate(reads):
            f.write(f'>{chrom}:{pos}|{i}\n')
            f.write(seq + '\n')

            reads_to_data[f'{chrom}:{pos}|{i}'] = {
                    'sequence': seq,
                    'sequence_base': base
                    }

            if i >= max_depth:
                break
    f.close()

    return reads_to_data
