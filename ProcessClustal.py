from Bio import AlignIO
import svgwrite as svg
from glob import glob
from os import path
from sys import argv

def process_seq_blocks(aln, set1, set2):
    set1_idxs = {i for i, rec in enumerate(aln) if rec.id.split('_')[0] in set1}
    set2_idxs = {i for i, rec in enumerate(aln) if rec.id.split('_')[0] in set2}
    current_block_type = ''
    current_block_start = 0
    current_block_stop = -1
    blocks = []
    for pos in range(len(aln[0])):
        bases_set1 = {aln[s, pos] for s in set1_idxs}
        bases_set2 = {aln[s, pos] for s in set2_idxs}

        if bases_set1 == bases_set2:
            next_block = 'IDENT'
        elif not bases_set1.intersection(bases_set2):
            next_block = 'DISJOINT'
        else:
            next_block = 'OTHER'

        if current_block_type == next_block:
            current_block_stop = pos + 1
        else:
            if current_block_type:
                blocks.append(
                    (current_block_start, current_block_stop, current_block_type)
                    )
            current_block_type = next_block
            current_block_start = pos
            current_block_stop = pos+1
    blocks.append((current_block_start, current_block_stop, current_block_type))
    return blocks



if __name__ == "__main__":
    pass

