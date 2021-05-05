#%% Libraries
import argparse

from pathlib import Path


import warnings
warnings.filterwarnings("ignore")

import pandas as pd
from pandas_plink import read_plink
import dask.array as da
import numpy as np

#%% Class definition
class FmapPrep():
    def __init__(self, method, save_dir,
                 ss_path,
                 G_dir, G_pattern,
                 ss_sep='\t',
                 ann_dir=None, ann_sep='\t', ann_pattern=None,
                 n_samples = None,
                 diag_eps=1e-4, off_eps=1e-5, debug=False):

        # Method
        self.method = method
        self.n_samples = n_samples
        method_list = ['riviera', 'paintor', 'finemap']
        assert self.method in method_list, f'method must be one of {method_list}'
        assert (method != 'finemap') or (isinstance(self.n_samples, int)), 'n_samples must be specified for finemap'

        # Save dir
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)

        # Sumstats
        self.ss_path = ss_path
        self.ss_sep = ss_sep

        # Individual-level data
        self.G_dir = Path(G_dir)
        self.G_pattern = G_pattern
        self.diag_eps = diag_eps
        self.off_eps = off_eps
        self.ld_ext = '.LD'

        # Annotations
        self.ann_dir = Path(ann_dir) if ann_dir is not None else None
        self.ann_sep = ann_sep
        self.ann_pattern = ann_pattern

        self.debug = debug


    def load_sumstats(self):
        sumstats = pd.read_csv(self.ss_path, sep=self.ss_sep)
        colnames = ['chr', 'snp', 'beta', 'se', 'locus']
        for colname in colnames:
            assert colname in sumstats.columns, f'{colname} not found in column names'
        return(sumstats)

    def format_sumstats(self, sumstats):
        sumstats['z'] = sumstats.beta/sumstats.se
        return(sumstats)


    def load_G(self, chr_n):
        file_path = str(self.G_dir/self.G_pattern.replace('*', str(chr_n)))
        bim, _, G = read_plink(file_path, verbose=False)
        G = G.T
        bim = bim[['snp', 'i']]
        return(bim, G)

    def load_ann(self, chr_n):
        '''
        expects: (P + 1) x (A + 1) table. Index should be rsids, header should be annotation name
        '''
        print('Loading annotations')
        file_path = str(self.ann_dir/self.ann_pattern.replace('*', str(chr_n)))
        ann = pd.read_csv(file_path, sep=self.ann_sep)
        return(ann)

    def format_ann(self, ann):
        ann = ann.drop(['CM', 'BP', 'CHR'], axis=1)
        ann = ann.rename({'SNP':'snp'}, axis=1)
        return(ann)

    def get_intersect_subset(self, *dfs):
        inter_snps = set.intersection(*[set(df.snp) for df in dfs])
        return([df[df.snp.isin(inter_snps)] for df in dfs])

    def alligned_locus_data(self, locus_ss, bim, ann):
        locus_snps = locus_ss.snp
        locus_bim = bim.set_index('snp').loc[locus_snps, :].reset_index()
        locus_ann = ann.set_index('snp').loc[locus_snps, :].reset_index()
        return(locus_ss, locus_bim, locus_ann)

    def calculate_LD(self, g):
        g = (g - g.mean(axis=0))/g.std(axis=0)
        LD = (da.dot(g.T, g)/ g.shape[0]).compute()
        return(LD)

    def correct_LD(self, LD):
        # Symmetry
        LD = (LD.T + LD)/2

        # Off diagonal elements
        LD = np.where(np.abs(LD) >= 1, np.sign(LD)*(1 - self.off_eps), LD)

        # Diagonal elements
        np.fill_diagonal(LD, 1 + self.diag_eps)
        return(LD)

    def save_LD(self, LD, filename):
        np.savetxt(self.save_dir/filename, LD, fmt='%f')

    def paintor_save(self, locus_ss, locus_ann, locus_save_name):
        locus_ss[['snp', 'z']].to_csv(self.save_dir/locus_save_name, index=None, sep=' ')
        locus_ann.drop('snp', axis=1).to_csv(self.save_dir/(locus_save_name + '.annotations'), index=None)
        return(locus_save_name)

    def finemap_save(self, locus_ss, locus_save_name):
        zname = locus_save_name + '.z'
        locus_z = locus_ss[['chr', 'snp', 'beta', 'se']]
        locus_z = locus_z.assign(allele1='A', allele2='G', position=9, maf=0.5)
        locus_z = locus_z.rename({'snp': 'rsid', 'chr': 'chromosome'}, axis=1)
        locus_z = locus_z[['rsid', 'chromosome', 'position',
                           'allele1', 'allele2', 'maf', 'beta', 'se']]
        locus_z.to_csv(self.save_dir/zname, index=None, sep=' ')
        return([zname, locus_save_name + self.ld_ext,
                *[locus_save_name + '.' + ext for ext in ('snp', 'config', 'cred', 'log') ]])

    def riviera_save(self, locus_ss, locus_ann, locus_save_name):
        self.paintor_save(locus_ss, locus_ann, locus_save_name)

    def paintor_tracker_save(self, tracker):
        tracker = pd.Series(tracker)
        tracker.to_csv(self.save_dir/'input.files', index=None, header=None)

    def finemap_tracker_save(self, tracker):
        tracker = pd.DataFrame(tracker,
                               columns=['z', 'ld', 'snp', 'config', 'cred', 'log'])
        tracker = tracker.assign(n_samples = self.n_samples)
        tracker.to_csv(self.save_dir/'master', index=None, sep=';')

    def prep(self):
        sumstats = self.load_sumstats()
        sumstats = self.format_sumstats(sumstats)

        tracker = []
        for i, (chr_n, chr_ss) in enumerate(sumstats.groupby('chr')):
            print(f'> Processing chromosome {chr_n}')
            bim, G = self.load_G(chr_n)

            if self.ann_dir is None or self.method == 'finemap':
                ann = pd.DataFrame({'snp':chr_ss.snp.unique(), 'base':1})
            else:
                ann = self.load_ann(chr_n)
                ann = self.format_ann(ann)

            chr_ss, bim, ann = self.get_intersect_subset(chr_ss, bim, ann)

            for locus_n, locus_ss in chr_ss.groupby('locus'):
                locus_save_name = f'Locus{int(locus_n)}'

                locus_ss, locus_bim, locus_ann = self.alligned_locus_data(locus_ss, bim, ann)


                g = G[:, locus_bim.i.to_list()]
                if not self.debug:
                    print(f'Calculating LD for locus {locus_n}')
                    LD = self.calculate_LD(g)
                    print('done')
                    LD = self.correct_LD(LD)
                else:
                    LD = np.eye(len(locus_ss))
                self.save_LD(LD, locus_save_name + self.ld_ext)

                if self.method == 'paintor':
                    tracker_info = self.paintor_save(locus_ss, locus_ann,
                                                     locus_save_name)
                    tracker.append(tracker_info)

                elif self.method == 'finemap':
                    tracker_info = self.finemap_save(locus_ss, locus_save_name)
                    tracker.append(tracker_info)

                elif self.method == 'riviera':
                    self.riviera_save(locus_ss, locus_ann, locus_save_name)

        if self.method == 'paintor':
            self.paintor_tracker_save(tracker)

        elif self.method == 'finemap':
            self.finemap_tracker_save(tracker)
        return(chr_ss)

#%% Main
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Utility script to prepare data for fine-mapping with different methods')

    # Positional arguments
    parser.add_argument('--method',
                        type=str,
                        help='name of fine-mapping method for which to prepare',
                        choices=['riviera', 'paintor', 'finemap'])

    parser.add_argument('--save_dir',
                        type=str,
                        help='path to directory to save final files')

    parser.add_argument('--ss_path',
                        type=str,
                        help='path to summary statistics file')

    parser.add_argument('--ss_sep',
                        type=str,
                        default='\t',
                        help='separator between summary statistics file columns')

    parser.add_argument('--G_dir',
                        type=str,
                        help='path to directory with 1000 Genome plink files (bim, bed, fam)')

    parser.add_argument('--G_pattern',
                        type=str,
                        help='pattern of 1000 Genome files with chromosome number replaced by "*". Ex: "1000G.EUR.QC.1.bim" -> "1000G.EUR.QC.*.bim"')

    parser.add_argument('--ann_dir',
                        type=str,
                        help='path to directory with annotation files')

    parser.add_argument('--ann_sep',
                        type=str,
                        default='\t',
                        help='separator between annotations file columns')

    parser.add_argument('--ann_pattern',
                        type=str,
                        help='pattern of annotation files with chromosome number replaced by "*". Ex: "baseline.1.annot.gz" -> "baseline.*.annot.gz"')

    parser.add_argument('--n_samples',
                        type=int,
                        help='Number of samples (individuals) used to calculate summary statistics (only necessary for FINEMAP)')

    parser.add_argument('--debug',
                        type=bool,
                        help='debug mode (identity LD)',
                        default=False,
                        choices=[True, False])

    args = parser.parse_args()

    #%%
    proc =  FmapPrep(method=args.method, save_dir=args.save_dir,
                     ss_path=args.ss_path, ss_sep=args.ss_sep,
                     G_dir=args.G_dir, G_pattern=args.G_pattern,
                     ann_dir=args.ann_dir, ann_sep=args.ann_sep, ann_pattern=args.ann_pattern,
                     debug=args.debug, n_samples=args.n_samples)

    proc.prep()
