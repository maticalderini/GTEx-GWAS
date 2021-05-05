#%% Libraries
from pathlib import Path
import pandas as pd
import numpy as np

import warnings
warnings.filterwarnings("ignore")

import argparse

#%% Class definition
class RefPrep():
    def __init__(self, snps_ref_path, gene_ref_path, save_dir):
        self.snps_ref_path = snps_ref_path
        self.gene_ref_path = gene_ref_path

        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)

    def load_snps_ref(self):
        print('Loading snps ref (may take a few minutes)')
        snps = pd.read_csv(self.snps_ref_path, '\t',
                           usecols=[0, 1, 2, 6], names=['chr', 'pos', 'id', 'rsid'],
                           low_memory=False)
        return(snps)

    def clean_snps(self, snps):
        print('Cleaning snps ref (may take a few minutes)')
        # chromosome column
        snps['chr'] = pd.to_numeric(snps['chr'], errors='coerce')
        snps = snps.dropna().astype({'chr':int})
        snps = snps.reset_index(drop=True)

        # Rsid present
        snps.loc[snps.rsid == '.', 'rsid'] = np.nan

        return(snps)

    def prep_snps(self):
        snps = self.load_snps_ref()
        snps = self.clean_snps(snps)
        return(snps)

    def load_gene_ref(self):
        print('Loading genes ref (may take a few minutes)')
        genes = pd.read_csv(self.gene_ref_path, header=None, sep='\t',
                             skiprows=6, usecols=[0, 2, 8], names=['chr', 'feature', 'meta'])
        return(genes)

    def clean_genes(self, genes):
        # clean chromosomes
        print('Cleaning genes ref (may take a few minutes)')
        genes['chr'] = pd.to_numeric(genes['chr'], errors='coerce')
        genes = genes.dropna()
        genes['chr'] = genes['chr'].astype(int)

        # Remove all non-gene
        genes = genes[genes.feature == 'gene'].drop('feature', axis=1)

        # Get meta info
        split_meta = genes.meta.str.split(';', expand=True)[[0, 4]]
        split_meta = split_meta.rename({0:'id', 4:'name'}, axis=1)
        split_meta = split_meta.apply(lambda col:col.str.extract('"([^"]*)"')[0])
        genes = pd.concat([genes.drop('meta', axis=1), split_meta], axis=1)

        genes = genes.reset_index(drop=True)
        return(genes)

    def prep_genes(self):
        genes = self.load_gene_ref()
        genes = self.clean_genes(genes)
        return(genes)

    def save_ref(self, ref, name):
        print(f'saving {name} (may take a few minutes)')
        ref.to_csv(self.save_dir/(name + '.txt.gz'),
                   sep='\t', index=None, compression='gzip')

    def save_snps_ref(self, ref_snps):
        self.save_ref(ref_snps, 'ref_snps')

    def save_genes_ref(self, ref_genes):
        self.save_ref(ref_genes, 'ref_genes')

    def prep_refs(self):
        snps = self.prep_snps()
        genes = self.prep_genes()

        self.save_snps_ref(snps)
        self.save_genes_ref(genes)
        return(snps, genes)


#%% Main
if __name__ == '__main__':

   parser = argparse.ArgumentParser(description='Utility script to prepare gtex reference files')

   # Positional arguments
   parser.add_argument('--gene_ref',
                       type=str,
                       help='raw gene reference file')

   parser.add_argument('--snps_ref',
                       type=str,
                       help='raw snp reference file')

   parser.add_argument('--save_dir',
                       type=str,
                       help='directory to save clean ref files')

   args = parser.parse_args()


   prep = RefPrep(snps_ref_path=args.snps_ref,
                  gene_ref_path=args.gene_ref,
                  save_dir=args.save_dir)

   prep.prep_refs()
