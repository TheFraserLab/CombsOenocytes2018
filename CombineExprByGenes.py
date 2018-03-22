import pandas as pd
from argparse import ArgumentParser, FileType
from sys import stdout
from tqdm import tqdm

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--outfile', '-o', default=stdout, type=FileType('w'))
    parser.add_argument('--use-gtf', default=False, action='store_true',
                        help='TRANSCRIPT_INFO is a GTF-like file')
    parser.add_argument('--use-gene-name', default=False, action='store_true',
                        help='Implies --use-gtf; Output is gene name instead of'
                        ' gene id')
    parser.add_argument('--cufflinks-tlst', default={})
    parser.add_argument('transcript_info', type=FileType('r'))
    parser.add_argument('expression_table')
    args =  parser.parse_args()
    if args.use_gene_name:
        args.use_gtf = True
    if args.cufflinks_tlst:
        args.cufflinks_tlst = dict([
            line.split()[:2] for line in open(args.cufflinks_tlst)
            ])
    return args


if __name__ == "__main__":
    args = parse_args()
    genes = set()
    tr_to_gene = {}

    if args.use_gtf:
        for line in tqdm(args.transcript_info):
            data = line.split('\t')
            annot = {e.split()[0]: e.split()[1].strip().strip('\'";')
                     for e in data[-1].strip().strip(';').split(';')}
            if args.use_gene_name:
                genes.add(annot['gene_name'])
                tr_to_gene[annot['transcript_id']] = annot['gene_name']
            else:
                genes.add(annot['gene_id'])
                tr_to_gene[annot['transcript_id']] = annot['gene_id']

    else:
        for line in args.transcript_info:
            gene, *transcripts = line.split()
            genes.add(gene)
            for transcript in transcripts:
                tr_to_gene[transcript] = gene

    expr_table = pd.read_table(args.expression_table, comment="#", index_col=0)

    out_table = pd.DataFrame(data=pd.np.nan, index=sorted(genes),
                             columns=expr_table.columns)

    pbar = tqdm(expr_table.iterrows(), total=len(expr_table))
    for tr, row in pbar:
        tr = args.cufflinks_tlst.get(str(tr), tr)
        gn = tr_to_gene[tr]
        if not row.count():
            continue
        if not out_table.loc[gn].count():
            out_table.loc[gn] = 0

        out_table.loc[gn] += row

    out_table.to_csv(args.outfile, sep='\t', )
    args.outfile.close()




