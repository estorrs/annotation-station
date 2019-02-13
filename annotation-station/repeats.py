import os
import re
from collections import defaultdict


CHROM_COLUMN = 5
START_COLUMN = 6
STOP_COLUMN = 7
NAME_COLUMN = 10
CLASS_COLUMN = 11
FAMILY_COLUMN = 12

def normalize_chrom(chrom):
    return re.sub(r'^chr', '', chrom)

def get_position_index(chrom, pos, bin_size=1000000):
    """Return index for the given genomic position"""
    # get rid of chr if present
    chrom = normalize_chrom(chrom)
    
    pos_tup = (chrom, int(int(pos) / bin_size))
    
    return pos_tup

class RepeatCollection(object):
    def __init__(self, bin_size=1000000):
        self.allias_to_tups = defaultdict(list)
        self.bin_size = bin_size

    def put_repeat(self, repeat_tup):
        """put repeat into collection

        repeat_tup - (chrom, start, stop, name, class, family)"""
        start_index = get_position_index(repeat_tup[0], repeat_tup[1],
                bin_size=self.bin_size)
        stop_index = get_position_index(repeat_tup[0], repeat_tup[2],
                bin_size=self.bin_size)

        self.allias_to_tups[start_index].append(repeat_tup)
        if start_index != stop_index:
            self.allias_to_tups[stop_index].append(repeat_tup)

    def get_repeat(self, chrom, pos):
        """Get repeat from collection"""
        position_index = get_position_index(chrom, pos,
                bin_size=self.bin_size)
        
        potentials = self.allias_to_tups[position_index]

        for repeat in potentials:
            is_in_range = int(pos) >= int(repeat[1]) and int(pos) <= int(repeat[2])
            is_same_chrom = normalize_chrom(repeat[0]) == normalize_chrom(chrom)
            if is_in_range and is_same_chrom:
                return repeat

        return None

def get_repeat_collection(repeat_table_fp):
    """Get repeat collection from file"""
    f = open(repeat_table_fp)

    # kill header
    f.readline()

    rc = RepeatCollection()
    for line in f:
        pieces = line.strip().split('\t')

        chrom = pieces[CHROM_COLUMN]
        start = pieces[START_COLUMN]
        stop = pieces[STOP_COLUMN]
        repeat_name = pieces[NAME_COLUMN]
        repeat_class = pieces[CLASS_COLUMN]
        repeat_family = pieces[FAMILY_COLUMN]

        rc.put_repeat((chrom, start, stop, repeat_name, repeat_class, repeat_family))

    return rc

# repeats table from ucsc table viewer
RP = get_repeat_collection(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/repeats_table.grch38.tsv'))


def get_repeat_by_position(chrom, pos):
    """Returns repeat if there is one at given position.

    returns (repeat_name, repeat_class, repeat_family)

    If no repeat is present at given position, then None is returned"""
    repeat = RP.get_repeat(chrom, pos)
    if repeat is not None:
        return repeat[3], repeat[4], repeat[5]
    return None
