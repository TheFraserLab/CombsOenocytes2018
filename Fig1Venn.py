import pandas as pd
import numpy as np
import re
from glob import glob
from tqdm import tqdm
from collections import defaultdict
from sys import argv
from os import path

true_index = lambda x: x.index[x]

def print_venn_stats(deseq_genes, f_genes, o_genes, fname):
    d = deseq_genes
    f = f_genes
    o = o_genes

    print("d", len(deseq_genes), "f", len(f_genes), "o", len(o_genes))
    print(
            "df", len(deseq_genes.intersection(f_genes)),
            "do", len(deseq_genes.intersection(o_genes)),
            "fo", len(o_genes.intersection(f_genes))
            )
    print("dfo", len(deseq_genes.intersection(f_genes).intersection(o_genes)))
    print('-'*30)
    print("A = ASE, B=sex, C=tissue")
    print('A !B !C', len(d.difference(f).difference(o)))
    print('B !A !C', len(f.difference(d).difference(o)))
    print('C !A !B', len(o.difference(f).difference(d)))

    print('A  B !C', len(d.intersection(f).difference(o)))
    print('A  C !B', len(d.intersection(o).difference(f)))
    print('B  C !A', len(f.intersection(o).difference(d)))

    print('A  B  C', len(d.intersection(f).intersection(o)))

    all_filters = d.intersection(f).intersection(o)
    tissue_spec_only = f.intersection(o).difference(d)
    ase_only = d.difference(all_filters)
    outf = open(fname, 'w')
    if len(all_filters):
        print(*all_filters, sep='\tallpass\n', end='\tallpass\n', file=outf)
    if len(tissue_spec_only):
        print(*tissue_spec_only, sep='\ttissue\n', end='\ttissue\n', file=outf)
    if len(ase_only):
        print(*ase_only, sep='\tase\n', end='\tase\n', file=outf)
    outf.close()

if __name__ == "__main__":
    fbtr = re.compile('FBtr[0-9]*')
    gn = re.compile('gene_name "([^"]*)"')

    dirname = argv[1]
    species = path.basename(dirname.strip('/'))

    orthologdb = {}
    orthologs = sorted(glob('prereqs/gene_orthologs*.tsv'))
    for fname in orthologs:
        with open(fname) as f:
            for line in f:
                if line.startswith('#'): continue
                if not line.strip(): continue
                data = line.split('\t')
                symbol = 'desatF' if data[1] == 'Fad2' else data[1]
                orthologdb[data[0]] = symbol
                orthologdb[data[6]] = symbol
                orthologdb[data[5]] = symbol

    fb_to_gn = {}
    for line in open('Reference/{}_good.gtf'.format(species)):
        the_fbtr = fbtr.findall(line)[0]
        the_gn = gn.findall(line)[0]
        fb_to_gn[the_fbtr] = the_gn

    deseq = pd.read_table(path.join(dirname, 'deseq_pvals.tsv'),
            index_col=0)

    for target_sex in ['male', 'female']:
        for target_tissue in ['oe', 'fb']:

            sleuth_sex = pd.read_table(path.join(
                dirname,
                'combined/{}_sleuth.tsv'.format(target_sex)),
                index_col=0 )
            sleuth_tissue = pd.read_table(
                    path.join(
                        dirname,
                        'combined/{}_sleuth.tsv'.format(target_tissue)
                        ),
                    index_col=0)



            tissue_obs_norm = pd.read_table(path.join(dirname,
                'combined/sleuth_{}_obs_norm.tsv'.format(target_tissue)))
            sex_obs_norm = pd.read_table(path.join(dirname,
                'combined/sleuth_{}_obs_norm.tsv'.format(target_sex)))

            self_vs_othersex = defaultdict(float)
            othersex = defaultdict(float)
            self_vs_othertissue = defaultdict(float)
            othertissue = defaultdict(float)

            for i in tqdm(tissue_obs_norm.index):
                row = tissue_obs_norm.loc[i]
                gene = row.target_id
                samp = row['sample']
                tpm = row.tpm
                ddl = [self_vs_othersex, othersex][not samp[2:].startswith(target_sex)]
                ddl[gene]+= (tpm)

            for i in tqdm(sex_obs_norm.index):
                row = sex_obs_norm.loc[i]
                gene = row.target_id
                samp = row['sample']
                tpm = row.tpm
                ddl = [self_vs_othertissue, othertissue][not samp[:2] == target_tissue]
                ddl[gene]+= (tpm)

            sleuth_sex['selfhigher'] = False
            sleuth_sex['self'] = 0
            sleuth_sex['other'] = 0
            sleuth_tissue['selfhigher'] = False
            sleuth_tissue['self'] = 0
            sleuth_tissue['other'] = 0

            for gene in tqdm(sleuth_sex.index):
                res = sex_obs_norm.query('target_id == "{}"'.format(gene))
                sleuth_sex.ix[gene, 'self'] = self_vs_othertissue[gene]
                sleuth_sex.ix[gene, 'other'] = othertissue[gene]
                sleuth_sex.ix[gene, 'selfhigher'] = (
                        (self_vs_othertissue[gene]) > (othertissue[gene])
                        )

            for gene in tqdm(sleuth_tissue.index):
                res = tissue_obs_norm.query('target_id == "{}"'.format(gene))
                sleuth_tissue.ix[gene, 'selfhigher'] = (
                        (self_vs_othersex[gene]) > (othersex[gene])
                        )

            sleuth_sex_gn = sleuth_sex.rename(index=fb_to_gn)
            sleuth_tissue_gn = sleuth_tissue.rename(index=fb_to_gn)

            cutoff = 1e-3
            deseq_genes = true_index(abs(deseq[target_tissue+target_sex]) > -np.log10(cutoff))
            f_genes = set(true_index((sleuth_sex_gn.qval < cutoff) &
                (sleuth_sex_gn.selfhigher)))
            o_genes = set(true_index((sleuth_tissue_gn.qval < cutoff) &
                (sleuth_tissue_gn.selfhigher)))

            print('-'*30)
            print(target_tissue, target_sex)
            print_venn_stats(deseq_genes, f_genes, o_genes,
                    fname=path.join(dirname, 'combined',
                        '{}{}_spec_genes.txt'.format(target_tissue, target_sex)
                        ))

