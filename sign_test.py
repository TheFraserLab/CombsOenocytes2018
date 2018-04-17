# Performs sign test using a .gmt file and an ase file (containing gene names and ase values)

from __future__ import print_function
import scipy.stats as stats
from argparse import ArgumentParser
import heapq as hq
import os

parser = ArgumentParser()
parser.add_argument('--min-genes', default=20, type=int)
parser.add_argument('--test-only', default=False,
                    help="Test only go terms with a given code in second column")
parser.add_argument('--drop-genes', default=None)
parser.add_argument('--print-header', default=False, action='store_true')
parser.add_argument('--print-counts', default=False, action='store_true')
parser.add_argument('--pseudocounts', default=1, type=int)
parser.add_argument('--ase-col', default=-1, type=int)
parser.add_argument('--ase-cutoff', default=0, type=float)
parser.add_argument('--translate-orthologs', default=False, action='store_true')
parser.add_argument('ase')
parser.add_argument('categories')
parser.add_argument('outfile')
parser.add_argument('cutoff', type=int)

args = parser.parse_args()
ase = args.ase
categories = args.categories
outfile = args.outfile
CUTOFF = args.cutoff

orthologdb = {}
if args.translate_orthologs:
    from glob import glob
    orthologdb = {}
    orthologs = sorted(glob(
        'prereqs/gene_orthologs*.tsv'))
    print(orthologs)
    for fname in orthologs:
        with open(fname) as f:
            for line in f:
                if line.startswith('#'): continue
                if not line.strip(): continue
                data = line.split('\t')
                symbol = data[1]
                orthologdb[data[0]] = symbol
                orthologdb[data[6]] = symbol
                orthologdb[data[5]] = symbol


top_genes = []
top_values = []
bottom_genes = []
bottom_values = []
top_heap = []
bottom_heap = []
heaps = [bottom_heap, top_heap]

drop_genes = set()
if args.drop_genes:
    for gene in open(args.drop_genes):
        drop_genes.add(gene.strip())

with open(ase, 'r') as ase:
        ase.readline()

        for line in ase:
                line = line.strip().split('\t')
                gene = line[0]
                gene = orthologdb.get(gene, gene)
                if gene in drop_genes: continue
                ase_value = line[args.ase_col]
                if ase_value == "NA":
                        continue
                else:
                        ase_value = float(ase_value)

                heap = heaps[ase_value > 0]
                if len(heap) < CUTOFF:
                    hq.heappush(heap, (abs(ase_value), gene))
                else:
                    hq.heappushpop(heap, (abs(ase_value), gene))

top_values, top_genes = zip(*top_heap)
bottom_values, bottom_genes = zip(*bottom_heap)


out = open(outfile, 'w')
if args.print_header:
    if args.print_counts:
        print("category\toddsratio\tpval\ttop_yes\ttop_no\tbottom_yes\tbottom_no\tty_genes\tby_genes", file=out)
    else:
        print("category\toddsratio\tpval", file=out)

tests = 0

with open(categories, 'r') as categories:
        for line in categories:
                line = line.strip().split('\t')
                category = line[0]
                genes = set(line[2:])
                if args.test_only and line[1] != args.test_only: continue

                bottom_yes = args.pseudocounts
                bottom_no = args.pseudocounts
                top_yes = args.pseudocounts
                top_no = args.pseudocounts

                ty = []
                by = []
                for i in top_genes:
                        if i in genes:
                                top_yes += 1
                                ty.append(i)
                        else:
                                top_no += 1
                for i in bottom_genes:
                        if i in genes:
                                bottom_yes += 1
                                by.append(i)
                        else:
                                bottom_no += 1

                if bottom_yes + top_yes >= args.min_genes:
                        tests += 1
                        oddsratio, pvalue = stats.fisher_exact([[bottom_yes, top_yes], [bottom_no, top_no]])
                        outlist = [category, str(oddsratio), str(pvalue)]
                        if args.print_counts:
                            outlist.extend([top_yes, top_no, bottom_yes,
                                            bottom_no, ','.join(ty),
                                            ','.join(by)])
                        print(*outlist, sep='\t', file=out)

out.close()

if tests > 0:
    print("Num top genes:", len(top_genes))
    print("Num bottom genes:", len(bottom_genes))
    print("Total number of tests: " + str(tests))
    print("Adjusted significance level: " + str(0.05 / tests))
else:
        print("Not enough genes tested.")
