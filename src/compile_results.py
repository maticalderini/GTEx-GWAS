#%% Libraries
from pathlib import Path
import pandas as pd

import argparse

#%% 
class GTExCompiler():
    def __init__(self, tissues_dir, save_dir):
        self.tissues_dir = Path(tissues_dir)
        
        self.save_dir = Path(save_dir)
        self.save_dir.parent.mkdir(parents=True, exist_ok=True)
    
    def get_group_ss(self, gene_group_dir):
        gene_group = gene_group_dir.name
        group_ss = pd.read_csv(gene_group_dir/(gene_group + '.ss.txt'), sep='\t')
        return(group_ss)
        
    
    def get_locus_results(self, locus_path):
        locus_file = locus_path.name
        gene_group_dir = locus_path.parent
        
        locus_n = int(locus_file[5:])
        locus_results = pd.read_csv(gene_group_dir/'results'/f'Locus{locus_n}.results', sep=' ')
        locus_results = locus_results.assign(locus=locus_n)
        return(locus_results)
        
    def merge_ss_results(self, group_ss, locus_results):
        '''
        Merging with left join, because some snps
        from group_ss where not fine-mapped (not part of a locus)
        '''
        merged = group_ss.merge(locus_results[['snp', 'locus', 'Posterior_Prob']], how='left',
                                left_on=['snp', 'locus'], right_on=['snp', 'locus'])
        merged = merged.fillna(0)
        return(merged)
        
    def concat_tissue(self, tissue_ss_list, tissue_name):
        tissue_ss = pd.concat(tissue_ss_list)
        tissue_ss = tissue_ss.assign(tissue=tissue_name).drop('locus', axis=1)
        return(tissue_ss)
        
    def save_tissue_results(self, tissue_results, tissue_name):
        tissue_results.to_csv(self.save_dir/(tissue_name + '_results.txt'), sep='\t')      
        
    def compile(self):
        print('Compiling GTEx results')
        for tissue_dir in self.tissues_dir.iterdir():
            tissue_name = tissue_dir.name
            tissue_ss_list = []
            
            print(f'Compiling results for {tissue_name}')
            for gene_group_dir in tissue_dir.glob('gene_group_*'):
                group_ss = self.get_group_ss(gene_group_dir)
                
                for locus_path in gene_group_dir.glob('Locus*[!.annotations|.LD]'):
                    locus_results = self.get_locus_results(locus_path)
                    merged = self.merge_ss_results(group_ss, locus_results)
                    tissue_ss_list.append(merged)
                    
            tissue_results = self.concat_tissue(tissue_ss_list, tissue_name)
            self.save_tissue_results(tissue_results, tissue_name)
            
class GWASCompiler():
    def __init__(self, results_dir, gwas_path, save_path):
        
        self.results_dir = Path(results_dir)
        self.gwas_path = Path(gwas_path)
        
        self.save_path = Path(save_path)
        self.save_path.parent.mkdir(parents=True, exist_ok=True)
    
    def load_all_results(self):
        gwas_result_files =  self.results_dir.glob('Locus*.results')
        gwas_results = pd.concat([pd.read_csv(result_file, delim_whitespace=True) for result_file in gwas_result_files])
        gwas_results = gwas_results.drop('z', axis=1)
        return(gwas_results)
        
    def load_gwas(self):
        gwas = pd.read_csv(self.gwas_path, delim_whitespace=True)
        return(gwas)
        
    def merge_results(self, gwas, gwas_results):
        gwas_results = gwas.merge(gwas_results, how='left', left_on='snp', right_on='snp')
        gwas_results['Posterior_Prob'] = gwas_results['Posterior_Prob'].fillna(0)
        return(gwas_results)
    
    def save_results(self, gwas_results):
        gwas_results.to_csv(self.save_path, sep='\t')
    
    def compile(self):
        print('Compiling GWAS results')
        gwas_results = self.load_all_results()
        gwas = self.load_gwas()
        
        gwas_results = self.merge_results(gwas, gwas_results)      
        
        self.save_results(gwas_results)
        return(gwas_results)

class ScoreCompiler():
    def __init__(self, tissue_results_dir, gwas_results_path, save_path):
        self.tissue_results_dir = Path(tissue_results_dir)
        self.gwas_results_path = Path(gwas_results_path)
        
        self.save_path = Path(save_path)
        self.save_path.parent.mkdir(parents=True, exist_ok=True)
        
    def load_tissue_results(self):
        tissue_files = self.tissue_results_dir.glob('*[!gwas]_results.txt')
        tissues_results = pd.concat(pd.read_csv(tissue_file, delim_whitespace=True) for tissue_file in tissue_files)
        return(tissues_results)
        
    def load_gwas_results(self):
        gwas_results = pd.read_csv(self.gwas_results_path, delim_whitespace=True)
        return(gwas_results)
    
    def merge_results(self, tissues_results, gwas_results):
        '''
        Inner join here because score only needed for overlapping genes
        '''
        merged_results = tissues_results.merge(gwas_results[['snp', 'Posterior_Prob']],
                                               left_on='snp', right_on='snp',
                                               suffixes=('_tissue', '_gwas'))
        merged_results['prod'] = merged_results['Posterior_Prob_tissue']*merged_results['Posterior_Prob_gwas']
        return(merged_results)
           
    def calculate_scores(self, merged_results):
        scores = merged_results.groupby(['tissue', 'gene'])['prod'].sum()
        return(scores)
        
    def calculate_tissue_score(self, merged_results):
        tissue_score = merged_results.groupby(['tissue'])['prod'].sum()
        return(tissue_score)
        
    def save_scores(self, scores):
        scores = scores.rename('dot')
        scores.to_csv(self.save_path, sep='\t')
        
    def compile(self):
        print('Compiling Scores')
        tissues_results = self.load_tissue_results()
        gwas_results = self.load_gwas_results()
        
        merged_results = self.merge_results(tissues_results, gwas_results)      
        scores = self.calculate_scores(merged_results)
        self.save_scores(scores)
        return(scores)

#%% 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Utility script to compile gtex fine-mapping results')
    
    # Positional arguments
    parser.add_argument('--gwas_results_dir',
                        type=str,
                        help='path to gwas fine-mapping results directory')
    
    parser.add_argument('--gwas_path',
                        type=str,
                        help='path to GWAS summary statistics file (not results)')
    
    parser.add_argument('--gwas_save_path',
                        type=str,
                        help='path to save gwas results')
    
    parser.add_argument('--tissues_dir',
                        type=str,
                        help="path to directory containintg tissues' reulsts folders")
    
    parser.add_argument('--tissues_save_dir',
                        type=str,
                        help='directory to save tissues compiled results')
    
    parser.add_argument('--scores_save_path',
                        type=str,
                        help='directory to save gene and tissue level results')
    
    args = parser.parse_args()

    gwas_compiler = GWASCompiler(results_dir=args.gwas_results_dir,
                                 gwas_path=args.gwas_path,
                                 save_path=args.gwas_save_path)
    gwas_compiler.compile()
    
    
    gtex_compiler = GTExCompiler(tissues_dir=args.tissues_dir,
                                 save_dir=args.tissues_save_dir)
    gtex_compiler.compile()
    
    
    score_compiler = ScoreCompiler(tissue_results_dir=args.tissues_save_dir,
                                   gwas_results_path=args.gwas_save_path,
                                   save_path=args.scores_save_path)
    score_compiler.compile()
    

#%%
#    base_dir = Path(r'C:\Users\USer1\Documents\Consulting\Mcgill\chronic_pain\data\drg')
#    temp_dir= base_dir/'temp'
#    tissues_dir = temp_dir/'gtex_tissues'
#    
#    gwas_path = base_dir/'proc'/'gwas.txt'
#    gwas_results_dir = temp_dir/'gwas'/'results'
#    
#    save_dir = base_dir/'proc'
#    
#    gwas_save_path = save_dir/'gwas_results.txt'
#    scores_save_path = save_dir/'scores.txt'
# 
#        
#    
#    gwas_compiler = GWASCompiler(results_dir=gwas_results_dir,
#                                 gwas_path=gwas_path,
#                                 save_path=gwas_save_path)
#    gwas_compiler.compile()
#    
#    
#    gtex_compiler = GTExCompiler(tissues_dir=tissues_dir,
#                                 save_dir=save_dir)
#    gtex_compiler.compile()
#    
#    
#    score_compiler = ScoreCompiler(tissue_results_dir=save_dir,
#                                   gwas_results_path=gwas_save_path,
#                                   save_path=scores_save_path)
#    scores = score_compiler.compile()
#    
    
