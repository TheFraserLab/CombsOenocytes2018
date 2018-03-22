""" Calculate a polarized log10 pvalue based on a binomial

Assume both replicates are comparable, and can just add together the reference
and alternate counts.
"""
import pandas as pd
import numpy as np
from sys import stdout
from scipy.stats import binom_test
from glob import glob
from argparse import ArgumentParser
from collections import Counter
import re

def parse_args():
    "Parse command line arguments"
    parser = ArgumentParser()

    parser.add_argument('--outfile', '-o', default=stdout)
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    out_data = {}

    alt_finder = re.compile('[ {]1: ([0-9]*)')
    ref_finder = re.compile('-1: ([0-9]*)')

    for tissue in ['fbfemale', 'fbmale', 'oefemale', 'oemale']:
        files_in_tissue = [f for f in args.files if tissue in f]
        ref_counts = Counter()
        alt_counts = Counter()
        all_counts = Counter()
        tissue_ref_counts = 0
        tissue_alt_counts = 0
        for fname in files_in_tissue:
            header = next(open(fname))
            tissue_ref_counts += int(ref_finder.findall(header)[0])
            tissue_alt_counts += int(alt_finder.findall(header)[0])

            tab = pd.read_table(fname, index_col=0, skip_blank_lines=True,
                    comment='#')
            for ix in tab.index:
                ref_counts[ix] += tab.loc[ix, 'ref_counts']
                alt_counts[ix] += tab.loc[ix, 'alt_counts']
                all_counts[ix] += tab.loc[ix, 'ref_counts'] + tab.loc[ix,
                        'alt_counts']
        out_data[tissue] = {
                ix: (( np.log10(binom_test(
                    [ref_counts[ix], alt_counts[ix]],
                    p=(tissue_ref_counts/(tissue_ref_counts +
                        tissue_alt_counts))
                    ))
                     * (-1 if alt_counts[ix] > ref_counts[ix] else 1))
                    if all_counts[ix] > 10 else 0)
                for ix in tab.index}
        out_data[tissue]['__base__p__'] = (tissue_ref_counts /
                (tissue_ref_counts + tissue_alt_counts))
    pd.DataFrame(out_data).sort_index().to_csv(args.outfile, sep='\t')
