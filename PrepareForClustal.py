from Bio import SeqIO
from argparse import ArgumentParser
from os import path

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--outdir', '-o', default='tmp')
    parser.add_argument('fname', nargs='+')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    seqs = {path.basename(path.splitext(f)[0])
            : {rec.id: rec for rec in SeqIO.parse(f, 'fasta')}
            for f in args.fname}
    species = sorted(seqs)
    recids = sorted(seqs[species[0]])

    for i, recid in enumerate(recids):
        outfn = path.join(args.outdir, 'region_{}.fasta'.format(i))
        with open(outfn, 'w') as outfh:
            for s in species:
                rec = seqs[s][recid]
                rec.id = '{}_{}'.format(s, i)
                SeqIO.write(rec, outfh, 'fasta')




