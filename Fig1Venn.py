import pandas as pd
import numpy as np
import re
from glob import glob
from tqdm import tqdm
from collections import defaultdict

true_index = lambda x: x.index[x]


if __name__ == "__main__":

    sleuth_fem = pd.read_table('analysis/sim/combined/female_sleuth.tsv',
            index_col=1, )
    sleuth_oe = pd.read_table('analysis/sim/combined/oenocyte_sleuth.tsv',
            index_col=1, )

    deseq = pd.read_table('analysis/sim/combined/oefemale_deseq.tsv',
            index_col=0)
    deseq_genes = true_index(deseq.padj < .01)
    oe_obs_norm = pd.read_table('analysis/sim/combined/sleuth_oe_obs_norm.tsv')
    female_obs_norm = pd.read_table('analysis/sim/combined/sleuth_female_obs_norm.tsv')

    oef_vs_m = defaultdict(float)
    oem = defaultdict(float)
    oef_vs_fb = defaultdict(float)
    fbf = defaultdict(float)

    for i in tqdm(oe_obs_norm.index):
        row = oe_obs_norm.loc[i]
        gene = row.target_id
        samp = row['sample']
        tpm = row.tpm
        ddl = [oef_vs_m, oem]['female' not in samp]
        ddl[gene]+= (tpm)

    for i in tqdm(female_obs_norm.index):
        row = female_obs_norm.loc[i]
        gene = row.target_id
        samp = row['sample']
        tpm = row.tpm
        ddl = [oef_vs_fb, fbf]['fb' in samp]
        ddl[gene]+= (tpm)

    sleuth_fem['oehigher'] = False
    sleuth_oe['femhigher'] = False

    for gene in tqdm(sleuth_fem.index):
        res = female_obs_norm.query('target_id == "{}"'.format(gene))
        sleuth_fem.ix[gene, 'oehigher'] = (oef_vs_fb[gene]) > (fbf[gene])

    for gene in tqdm(sleuth_oe.index):
        res = oe_obs_norm.query('target_id == "{}"'.format(gene))
        sleuth_oe.ix[gene, 'femhigher'] = (oef_vs_m[gene]) > (oem[gene])


    fbtr = re.compile('FBtr[0-9]*')
    gn = re.compile('gene_name "([^"]*)"')

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
    for line in open('Reference/sim_good.gtf'):
        the_fbtr = fbtr.findall(line)[0]
        the_gn = gn.findall(line)[0]
        fb_to_gn[the_fbtr] = the_gn

    sleuth_f_gn = sleuth_fem.rename(index=fb_to_gn)
    sleuth_o_gn = sleuth_oe.rename(index=fb_to_gn)

    f_genes = set(true_index((sleuth_f_gn.qval < .01) & (sleuth_f_gn.oehigher)))
    o_genes = set(true_index((sleuth_o_gn.qval < .01) & (sleuth_o_gn.femhigher)))

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
    print('A !B !C', len(d.difference(f).difference(o)))
    print('B !A !C', len(f.difference(d).difference(o)))
    print('C !A !B', len(o.difference(f).difference(d)))

    print('A  B !C', len(d.intersection(f).difference(o)))
    print('A  C !B', len(d.intersection(o).difference(f)))
    print('B  C !A', len(f.intersection(o).difference(d)))

    print('A  B  C', len(d.intersection(f).intersection(o)))
