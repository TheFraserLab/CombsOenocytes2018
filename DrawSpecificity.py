import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from matplotlib import pyplot as mpl
from matplotlib.style import use
from scipy.stats import chi2
use('grayscale')

true_index = lambda x: x.index[x]

colors = {
        'allpass': (1,0,0, 1),
        'tissue': (0, 0, 1, 1),
        'ase': (0, 1.0, 0, 1),
        'other': (0, 0, 0, 0.2)
        }

startswith = lambda x: lambda y: y.startswith(x)
notstartswith = lambda x: lambda y: not y.startswith(x)

sw = startswith
nsw = notstartswith

arrowprops_below = dict(
        arrowstyle = "->",
        #connectionstyle = "arc3,angleA=0,angleB=-90,rad=3"
        connectionstyle = "arc3,rad=-0.3"
        )
arrowprops_above = arrowprops_below.copy()
arrowprops_above['connectionstyle'] = "arc3,angleA=0,angleB=90,rad=3"
arrowprops_above['connectionstyle'] = "arc3,rad=0.3"

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--file-types', '-t', nargs='+', default=['', '.png'])
    parser.add_argument('--try-orthologs', default=False)
    parser.add_argument('--pvals', default=None)
    parser.add_argument('--specificities')
    parser.add_argument('expr_file')
    parser.add_argument('ase_file')
    parser.add_argument('outdir')
    args = parser.parse_args()
    return args

def fishers_method(log10pvals):
    if len(log10pvals) == 1:
        return log10pvals[0]
    signs = set(np.sign(log10pvals))
    df = 2 * len(log10pvals)
    if len(signs) > 1:
        return 0
    elif signs == {-1}:
        chi2sum = sum(-log10pvals) * np.log(10)
        return np.log10(chi2.sf(chi2sum, df))
    elif signs == {1}:
        chi2sum = sum(log10pvals) * np.log(10)
        return -np.log10(chi2.sf(chi2sum, df))

if __name__ == "__main__":
    args = parse_args()

    expr = pd.read_table(args.expr_file, index_col=0, na_values=['---'])
    ase  = pd.read_table(args.ase_file,  index_col=0, na_values=['---'])

    if args.try_orthologs:
        from glob import glob
        orthologdb = {}
        orthologs = sorted(glob(
            os.path.join(args.try_orthologs, 'gene_orthologs*.tsv')))
        print(orthologs)
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
    else:
        orthologdb = {}


    if args.pvals is not None:
        args.pvals = pd.read_table(args.pvals, index_col=0, na_values=['---'])

    in_both = expr.index.intersection(ase.index)
    expr = expr.loc[in_both]
    ase  =  ase.loc[in_both]

    os.makedirs(args.outdir, exist_ok=True)

    for sex in ['male', 'female']:
        for tissue in ['oe', 'fb']:
            spec_file = os.path.join(args.specificities,
                    tissue+sex+'_spec_genes.txt')
            specificities = dict(line.strip().split()
                    for line in open(spec_file))
            specificity = (expr.T.select(sw(tissue+sex)).mean()
                           / expr.T.select(nsw(tissue+sex)).max().clip(0.1))
            tissue_ase = ase.T.select(sw(tissue+sex)).mean()

            mpl.figure(figsize=(4,3))
            if args.pvals is not None:
                asep = pd.Series({ix:
                    (args.pvals.loc[ix].select(sw(tissue+sex)).mean()
                    if ix in args.pvals.index
                    else 0)
                    for ix in tissue_ase.index
                    })
                #asep = asep.clip(-30, 30)
                asep = asep.clip(-350, 350)
                ix = asep.abs().sort_values().index
                asep = asep.loc[ix]
                specificity = specificity.loc[ix]
                tissue_ase = 50 * (tissue_ase.loc[ix] + 1)

                vmin = -(asep.abs().max())
                vmax = asep.abs().max()
            else:
                asep = pd.Series([(0, 0, 0, .1) for c in tissue_ase.index])
                vmin = vmax = 1
            cs = [ colors[specificities.get(gene, 'other')]
                    for gene in specificity.index]
            #assert tissue+sex != 'oefemale'
            mpl.scatter(specificity, asep, c=cs,
                    s=5,
                    edgecolors='none',
                    cmap=mpl.cm.Greys,
                    vmin=-1, vmax=1)
            mpl.xlabel('Specificity')
            mpl.ylabel('$-\log_{10} p $')
            low, hi = mpl.xlim()
            hi = max(hi, 1e3)
            ax = mpl.gca()
            ax.set_xscale('log', basex=10)
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position('top')
            mpl.xlim(0.095, hi*1.5)
            mpl.ylim(vmin*1.1, vmax*1.1)
            ytickvals = mpl.yticks()[0][1:-1]
            ytickvals = ytickvals.astype(int)
            mpl.yticks(ytickvals, np.abs(ytickvals),)

            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_bounds(0, 100)
            ax.spines['left'].set_bounds(-30, 30)
            ax.spines['left'].set_bounds(vmin, vmax)
            ax.spines['top'].set_bounds(0.1, hi)
            mpl.tight_layout()
            for ftype in args.file_types:
                mpl.savefig(os.path.join(args.outdir, tissue+sex+ftype))

            is_interesting = ((specificity > 3)
                    #& (abs(tissue_ase) > .2)
                    & (abs(asep) > 3))
            interesting_genes = tissue_ase.loc[is_interesting].sort_values().index
            interesting_genes = {gene for gene in specificities
                    if specificities[gene] == 'allpass'}
            ig = specificity.loc[is_interesting].sort_values().index
            interesting_genes = asep.loc[ig].sort_values(kind='mergesort').index
            interesting_genes = {gene for gene in specificities
                    if specificities[gene] == 'allpass'}
            print(sorted(interesting_genes),len(interesting_genes),  sep='\n')
            mpl.scatter(x=specificity.loc[interesting_genes],
                    y=asep.loc[interesting_genes],
                    s=5,
                    c='red',
                    edgecolors='none',
                    )
            for i, gene in zip(np.linspace(vmin, vmax, len(interesting_genes),
                                        endpoint=True),
                               interesting_genes):
                newgene = orthologdb.get(gene, gene)
                print(tissue, sex, gene, newgene, specificity.loc[gene],
                      tissue_ase.loc[gene], asep.loc[gene],
                      '{:0.2g}'.format(10**(-abs(args.pvals.loc[gene][tissue+sex]))),
                      sep='\t')
                #if newgene not in {'CG8534', 'bond', 'eloF', 'Fad2', }:
                    #continue
                if i > tissue_ase.loc[gene]:
                    arrowprops = arrowprops_above
                else:
                    arrowprops = arrowprops_below
                mpl.text(
                        specificity.loc[gene]*1.1,
                        asep.loc[gene],
                        newgene,
                        fontsize=8,
                        horizontalalignment='left',
                        verticalalignment='baseline' if asep.loc[gene] > 0 else 'top',
                        )
#                mpl.annotate(newgene,
#                             #[specificity.loc[gene], tissue_ase.loc[gene]],
#                             [specificity.loc[gene], asep.loc[gene]],
#                             [4e2, i],
#                             arrowprops=arrowprops,
#                             annotation_clip=False,
#                             horizontalalignment='left',
#                             fontsize=8,
#                             )
            #cbar.set_clim(*clim)
            for ftype in args.file_types:
                mpl.savefig(
                        os.path.join(args.outdir, tissue+sex+'_annot'+ftype),
                        dpi=900,
                        )
            mpl.close()
            if args.pvals is not None:
                mpl.figure()
                mpl.hist(asep.clip(-30, 30), bins=np.arange(-30,30.5,0.5),
                        range=(-30,30))
                mpl.savefig(os.path.join(args.outdir,
                    'hist'+tissue+sex+args.file_types[0]))
