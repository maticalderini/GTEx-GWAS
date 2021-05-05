from pathlib import Path

import dask.dataframe as dd
import pandas as pd
import numpy as np
import argparse

#%%

class GTExPreProc():
    def __init__(self, save_dir, gene_ref_path, snps_ref_path,
                 gtex_path, gwas_path,
                 gtex_sep='\t', gwas_sep='\t', save_sep='\t',
                 ngroups=100, sig_th=5e-6):
        # Ref Paths
        self.gene_ref_path = gene_ref_path
        self.snps_ref_path = snps_ref_path
        
        #Data paths
        self.gtex_path = gtex_path
        self.gtex_sep = gtex_sep
        
        self.gwas_path = gwas_path
        self.gwas_sep = gwas_sep
        
        # Save dir
        self.save_sep = save_sep
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)
        
        # Other
        self.ngroups = ngroups
        self.sig_th = sig_th        

    
    def load_refs(self):
        print('Loading ref files (may take a few minutes)')
        ref_snps = pd.read_csv(self.snps_ref_path, sep='\t').dropna()
        ref_genes = pd.read_csv(self.gene_ref_path, sep='\t')
        return(ref_snps, ref_genes)
    
    def load_gtex(self):
        gtex = dd.read_csv(self.gtex_path, sep=self.gtex_sep, usecols=[0, 1, 6])
        return(gtex)
    
    def load_gwas(self):
        gwas = pd.read_csv(self.gwas_path, sep=self.gwas_sep)
        return(gwas)
    
    def sub_overlap_genes(self, ref_snps, gwas, gtex):
        snps_overlap = ref_snps[ref_snps.rsid.isin(gwas.snp)]
        gtex = gtex.merge(snps_overlap[['chr', 'id']],
                          left_on='variant_id', right_on='id', how='inner')

        gtex = gtex.drop(['variant_id', 'id'], axis=1)
        return(gtex)
        
    def sub_sig_genes(self, gtex):
        is_sig_gene = gtex.groupby('gene_id')['pval_nominal'].min() < self.sig_th
        sig_genes = is_sig_gene[is_sig_gene].index.to_list()
        gtex = gtex[gtex.gene_id.isin(sig_genes)]
        
        gtex = gtex.drop_duplicates('gene_id')[['chr', 'gene_id']].sort_values('chr').reset_index(drop=True)
        return(gtex)
    
    def add_gene_name(self, gtex, ref_genes):
        gtex = gtex.merge(ref_genes[['id', 'name']],
                          left_on='gene_id', right_on='id', how='inner')
        gtex = gtex.drop('gene_id', axis=1)
        return(gtex)    
    
    def split_gene_groups(self, genes):
        gene_groups = np.array_split(genes, min(len(genes), self.ngroups))
        return(gene_groups)
    
    def save_gene_groups(self, gene_groups):
        print('Saving gene groups')
        for i, gene_group in enumerate(gene_groups):
            gene_group.to_csv(self.save_dir/f'gene_group_{i}.txt',
                              sep=self.save_sep, index=None)
        return
    
    def prepare_gtex(self):
        ref_snps, ref_genes = self.load_refs()
        gwas = self.load_gwas()
        gtex = self.load_gtex()
        
        gtex = self.sub_overlap_genes(ref_snps, gwas, gtex)
        
        print('Computing GTEx-GWAS overlapping genes')
        gtex = gtex.compute()
        
        gtex = self.sub_sig_genes(gtex)
       
        gtex = self.add_gene_name(gtex, ref_genes)

        gene_groups = self.split_gene_groups(gtex)
        self.save_gene_groups(gene_groups)
        return(gene_groups)

if __name__ == '__main__':
        
#    parser = argparse.ArgumentParser(description='Utility script to prepare gtex reference files')
#    
#    # Positional arguments
#    parser.add_argument('--save_dir',
#                        type=str,
#                        help='directory to save gene groups')
#    
#    parser.add_argument('--gene_ref_path',
#                        type=str,
#                        help='path to gene reference file')
#    
#    parser.add_argument('--snps_ref_path',
#                        type=str,
#                        help='path to snps reference file')
#    
#    parser.add_argument('--gtex_path',
#                        type=str,
#                        help='path to tissue gtex data')
#    
#    parser.add_argument('--gwas_path',
#                        type=str,
#                        help='path to gwas data')
#    
#    parser.add_argument('--ngroups',
#                        type=int,
#                        default=100,
#                        help='number of gene groups')
#        
#    parser.add_argument('--save_sep',
#                        type=str,
#                        default='\t',
#                        help='gene group column separator')
#    
#    args = parser.parse_args()
#    
#    
#    prep = GTExPreProc(save_dir=args.save_dir,
#                       gene_ref_path=args.gene_ref_path,
#                       snps_ref_path=args.snps_ref_path,
#                       gtex_path=args.gtex_path,
#                       gwas_path=args.gwas_path,
#                       ngroups=args.ngroups,
#                       save_sep=args.save_sep)
#    
#    prep.prepare_gtex()
   
    
#%%    
    
    # Parameters
    data_dir = Path(r'C:\Users\USer1\Documents\Consulting\Mcgill\chronic_pain\data\drg')
    ref_dir = data_dir/'proc'/'GTEx_ref'
    
    gene_ref_path = ref_dir/'ref_genes.txt.gz'
    snps_ref_path = ref_dir/'ref_snps.txt.gz'
    
    gtex_path = data_dir/'raw'/'Brain_Nucleus_accumbens_basal_ganglia.allpairs.txt'
    gtex_sep = '\t'
    
    gwas_path = data_dir/'proc'/'gwas.txt'
    
    ngroups=25
    
    save_dir = data_dir/'temp'/'gtex_tissues'/'BNABG'
    save_sep = '\t'
    
    
    
    prep = GTExPreProc(save_dir=save_dir, gene_ref_path=gene_ref_path, snps_ref_path=snps_ref_path,
                       gtex_path=gtex_path, gwas_path=gwas_path,
                       ngroups=ngroups, save_sep=save_sep)
    
    gtex = prep.prepare_gtex()
    
    # ref_snps, ref_genes = prep.load_refs()
    # gwas = prep.load_gwas()
    # gtex = prep.load_gtex()
    
    # GTEx ref
    # prep = RefPrep(snps_ref_path=snps_ref_path, gene_ref_path=gene_ref_path)
    
    # snps = prep.prep_snps()  
    # genes = prep.prep_genes()
    
    # # GWAS
    # gwas = pd.read_csv(gwas_path, '\t')
    
    # # Gtex Prep
    # print('Loading GTEx')
    # gtex = dd.read_csv(gtex_path, sep=gtex_sep)


    # #
    # snps_overlap = snps[snps.rsid.isin(gwas.snp)]
    # gtex_overlap = gtex[['gene_id', 'variant_id']].merge(snps_overlap[['chr', 'id']],
    #                                                      left_on='variant_id', right_on='id', how='inner')
    # gtex_overlap = gtex_overlap.drop_duplicates('gene_id')
    # gtex_overlap = gtex_overlap.drop(['variant_id', 'id'], axis=1)
    
    # gtex_overlap = gtex_overlap.merge(genes[['id', 'name']], left_on='gene_id', right_on='id', how='inner')
    # gtex_overlap = gtex_overlap.drop('gene_id', axis=1)
    
    # # .loc[gtex.variant_id.isin(snps_overlap.id), ['gene_id', 'variant_id']]
    # print('Computing overlap (may take a few minutes)')
    # gtex_overlap_c = gtex_overlap.compute()
