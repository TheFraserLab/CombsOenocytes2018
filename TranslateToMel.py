from argparse import ArgumentParser, FileType
from glob import glob
import re

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--delimiter', '-d', default='\t')
    parser.add_argument('gtf_file', type=FileType('r'))
    parser.add_argument('in_file', type=FileType('r'))
    parser.add_argument('out_file', type=FileType('w'))
    return parser.parse_args()


if __name__ == "__main__":
    fbtr = re.compile('FBtr[0-9]*')
    gn = re.compile('gene_name "([^"]*)"')

    args = parse_args()
    orthologdb = {}
    orthologs = sorted(glob('prereqs/gene_orthologs*.tsv'))
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

    fb_to_gn = {}
    for line in args.gtf_file:
        the_fbtr = fbtr.findall(line)[0]
        the_gn = gn.findall(line)[0]
        fb_to_gn[the_fbtr] = the_gn

    for line in args.in_file:
        data = line.strip().split('\t')
        for i in range(len(data)):
            gn = fb_to_gn.get(data[i], data[i])
            data[i] = orthologdb.get(gn, gn)
        print(*data, sep='\t', file=args.out_file)
    args.out_file.close()

