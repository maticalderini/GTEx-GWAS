from pathlib import Path

import pandas as pd
import dask.dataframe as dd

import numpy as np

from utils import z2p

import argparse

#%%
class GTExProc():
    def __init__(self, gene_group_path, gtex_path, snps_ref_path, save_path,
                 gene_group_sep='\t', gtex_sep='\t', save_sep='\t',
                 sig_th=5e-6, delta=50e3):

        #Gene group
        self.gene_group_path = gene_group_path
        self.gene_group_sep = gene_group_sep

        #Gtex
        self.gtex_path = gtex_path
        self.gtex_sep = gtex_sep

        #Snps ref
        self.snps_ref_path = snps_ref_path

        #Save path
        self.save_sep = save_sep
        self.save_path = Path(save_path)
        self.save_path.parent.mkdir(parents=True, exist_ok=True)

        #Other
        self.sig_th = sig_th
        self.delta = delta


    def load_gene_group(self):
        gene_group = pd.read_csv(self.gene_group_path, sep=self.gene_group_sep)
        return(gene_group)

    def load_gtex(self):
        gtex = dd.read_csv(self.gtex_path, sep=self.gtex_sep)
        return(gtex)

    def format_gtex(self, gtex):
        gtex_cols = {'gene_id':'gene', 'variant_id':'snp_id',
             'pval_nominal':'p', 'slope':'beta', 'slope_se':'se'}
        gtex = gtex[list(gtex_cols.keys())]
        gtex.columns = list(gtex_cols.values())
        return(gtex)

    def get_group_gtex(self, gene_group, gtex):
        gtex = gtex.merge(gene_group, left_on='gene', right_on='id', how='inner')
        gtex = gtex.drop(['gene', 'id'], axis=1).rename(columns={'name':'gene'})
        gtex = gtex[['chr', 'gene', 'snp_id', 'beta', 'se']]
        print('computing group GTEx (may take a few minutes)')
        gtex = gtex.compute()
        return(gtex)

    def load_snps_ref(self):
        print('Loading snps ref (may take a few minutes)')
        snps_ref = pd.read_csv(self.snps_ref_path, sep='\t')
        return(snps_ref)

    def get_rsids(self, gtex, snps_ref):
        gtex = gtex.merge(snps_ref[['id', 'pos', 'rsid']],
                          left_on='snp_id', right_on='id', how='inner')
        gtex = gtex.drop(['id', 'snp_id'], axis=1).rename(columns={'rsid':'snp'})
        return(gtex)

    def sort_gtex(self, gtex):
        gtex = gtex.sort_values(['chr', 'gene', 'pos'])
        gtex = gtex[['chr', 'gene', 'snp', 'pos', 'beta', 'se']]
        return(gtex)

    def calc_log10p(self, gtex):
        z = gtex.beta/gtex.se
        p = z2p(z)
        log10p = -np.log10(p)
        return(log10p)

    def merge_loci(self, gene_stats):
        loci = gene_stats.locus.unique()
        for i in range(1, len(loci)):
            if (gene_stats.loc[gene_stats.locus==loci[i], 'pos'].min() - gene_stats.loc[gene_stats.locus==loci[i-1], 'pos'].max()) <= self.delta:
                gene_stats.loc[gene_stats.locus==loci[i-1], 'locus'] = loci[i]
        gene_stats.loc[:, 'locus'] = gene_stats.locus.astype('category').cat.codes + 1
        return(gene_stats)

    def create_sig_window_loci(self, gene_stats):
        gene_stats = gene_stats.assign(locus=np.nan)
        sig_pos = gene_stats.loc[gene_stats.log10p > -np.log10(self.sig_th), 'pos']
        for i, pos in enumerate(sig_pos, 1):
            gene_stats.loc[np.abs(gene_stats.pos - pos) <= self.delta, 'locus'] = i
        gene_stats = gene_stats.dropna().reset_index(drop=True)

        gene_stats = self.merge_loci(gene_stats)
        return(gene_stats)

    def correct_cross_gene_loci(self, gtex):
        new_loci = gtex.groupby(['gene', 'locus'])['chr'].unique().reset_index().rename_axis('new_locus').reset_index().drop('chr', axis=1)
        new_loci.index = range(1, len(new_loci) + 1)
        gtex = gtex.merge(new_loci, left_on=['gene', 'locus'], right_on=['gene', 'locus'])
        gtex = gtex.drop('locus', axis=1).rename(columns={'new_locus':'locus'})
        return(gtex)

    def divide_into_loci(self, gtex):
        print('Dividing into loci')
        gtex['log10p'] = self.calc_log10p(gtex)

        gtex = pd.concat([self.create_sig_window_loci(gene_stats) for gene, gene_stats in gtex.groupby('gene')])
        gtex = gtex.reset_index(drop=True)
        gtex = self.correct_cross_gene_loci(gtex)

        gtex = gtex.drop('log10p', axis=1)
        return(gtex)

    def save_sumstats(self, gtex):
        gtex.to_csv(self.save_path, sep=self.save_sep, index=None)

    def proc_gtex(self):
        gene_group = self.load_gene_group()

        gtex = self.load_gtex()
        gtex = self.format_gtex(gtex)

        gtex = self.get_group_gtex(gene_group, gtex)

        snps_ref = self.load_snps_ref()

        gtex = self.get_rsids(gtex, snps_ref)
        gtex = self.sort_gtex(gtex)

        gtex = self.divide_into_loci(gtex)

        self.save_sumstats(gtex)
        return(gtex)

#%%
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Utility script to prepare gtex reference files')

    # Positional arguments
    parser.add_argument('--gene_group_path',
                        type=str,
                        help='directory to gene group file')

    parser.add_argument('--snps_ref_path',
                        type=str,
                        help='path to snps reference file')

    parser.add_argument('--gtex_path',
                        type=str,
                        help='path to tissue gtex data')

    parser.add_argument('--save_path',
                        type=str,
                        help='path to save summary statistics')

    args = parser.parse_args()


    proc = GTExProc(gene_group_path=args.gene_group_path,
                    gtex_path=args.gtex_path,
                    snps_ref_path=args.snps_ref_path,
                    save_path=args.save_path)

    gtex = proc.proc_gtex()
