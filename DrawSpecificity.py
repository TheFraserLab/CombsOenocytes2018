import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from matplotlib import pyplot as mpl
from matplotlib.style import use
from scipy.stats import chi2
use('grayscale')



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
                    symbol = data[1]
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
            specificity = (expr.T.select(sw(tissue+sex)).mean()
                           / expr.T.select(nsw(tissue+sex)).max().clip(0.1))
            tissue_ase = ase.T.select(sw(tissue+sex)).mean()

            mpl.figure()
            if args.pvals is not None:
                color = pd.Series({ix: (fishers_method(args.pvals.loc[ix].select(sw(tissue+sex)))
                    if ix in args.pvals.index
                    else 0)
                    for ix in tissue_ase.index
                    }).clip(-30, 30)
                ix = color.abs().sort_values().index
                color = color.loc[ix]
                specificity = specificity.loc[ix]
                tissue_ase = 50 * (tissue_ase.loc[ix] + 1)

                vmin = -(color.abs().max())
                vmax = color.abs().max()
            else:
                color = pd.Series([(0, 0, 0, .1) for c in tissue_ase.index])
                vmin = vmax = 1
            mpl.scatter(specificity, tissue_ase, c=color, cmap=mpl.cm.coolwarm,
                    vmin=vmin, vmax=vmax)
            mpl.xlabel('Specificity')
            mpl.ylabel('% D. sechellia')
            low, hi = mpl.xlim()
            hi = max(hi, 1e3)
            ax = mpl.gca()
            cbar = mpl.colorbar(orientation='horizontal')
            clim = cbar.get_clim()
            ax.set_xscale('log', basex=10)
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position('top')
            #ax.xaxis.set_tick_params(which='both', labeltop='on', labelbottom='off', top='on',
                    #bottom='off')
            mpl.xlim(0.095, hi*1.5)
            mpl.ylim(-10, 120)
            mpl.yticks(np.arange(0, 101, 20))
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_bounds(0, 100)
            ax.spines['top'].set_bounds(0.1, hi)
            for ftype in args.file_types:
                mpl.savefig(os.path.join(args.outdir, tissue+sex+ftype))

            is_interesting = ((specificity > 3)
                    #& (abs(tissue_ase) > .2)
                    & (abs(color) > 3))
            interesting_genes = tissue_ase.loc[is_interesting].sort_values().index
            for i, gene in zip(np.linspace(0, 110, len(interesting_genes),
                                        endpoint=True),
                               interesting_genes):
                newgene = orthologdb.get(gene, gene)
                print(tissue, sex, gene, newgene, specificity.loc[gene],
                      tissue_ase.loc[gene], sep='\t')
                if i > tissue_ase.loc[gene]:
                    arrowprops = arrowprops_above
                else:
                    arrowprops = arrowprops_below
                mpl.annotate(newgene,
                             [specificity.loc[gene], tissue_ase.loc[gene]],
                             [4e2, i],
                             arrowprops=arrowprops,
                             annotation_clip=False)
            cbar.set_clim(*clim)
            for ftype in args.file_types:
                mpl.savefig(os.path.join(args.outdir, tissue+sex+'_annot'+ftype))
            mpl.close()
            if args.pvals is not None:
                mpl.figure()
                mpl.hist(color.clip(-30, 30), bins=np.arange(-30,30.5,0.5),
                        range=(-30,30))
                mpl.savefig(os.path.join(args.outdir,
                    'hist'+tissue+sex+args.file_types[0]))
