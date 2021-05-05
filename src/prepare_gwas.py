#%% Libraries
from pathlib import Path

import pandas as pd
import warnings
warnings.filterwarnings("ignore")

import numpy as np
from utils import z2p

import argparse

#%% Function definitions
class GWASPrep():
    def  __init__(self, ss_path, breaks_path, save_path,
                  ss_sep='\t', breaks_sep='\t', delta=0.5e6):
        '''
        delta: float window-size in base pairs (unlike command line which accepts it in Mb)
        '''
        # Sumstats
        self.ss_path = Path(ss_path)
        self.ss_sep = ss_sep
        
        #Breaks
        self.breaks_path = Path(breaks_path)
        self.breaks_sep = breaks_sep
        
        # Save Path
        self.save_path = Path(save_path)
        self.save_path.parent.mkdir(parents=True, exist_ok=True)
        
        #Other
        self.delta = delta
    
    def clean_sumstats(self, sumstats):
        sumstats = sumstats.dropna()
        sumstats = sumstats.sort_values('p_value').drop_duplicates(subset='SNP', keep='first')
        sumstats = sumstats.sort_values(['CHR', 'BP'])
        return(sumstats)
        
    def format_sumstats(self, sumstats):
        sumstats.columns = map(str.lower, sumstats.columns)
        sumstats = sumstats[['snp', 'chr', 'beta', 'se', 'bp']]
        
        sumstats['z'] = sumstats.beta/sumstats.se
        sumstats['p'] = z2p(sumstats.z)
        sumstats['log10p'] = -np.log10(sumstats.p)
        return(sumstats)
    
    def load_sumstats(self):
        print('Loading sumstats')
        sumstats = pd.read_csv(self.ss_path, sep=self.ss_sep)
        
        print('Cleaning and formating sumstats')
        sumstats = self.clean_sumstats(sumstats)
        sumstats = self.format_sumstats(sumstats)
        return(sumstats)
        
    
    def load_breaks(self):
        breaks = pd.read_csv(self.breaks_path, sep=self.breaks_sep, usecols=['GenomicLocus', 'chr', 'start', 'end']).set_index('GenomicLocus')
        return(breaks)
    
    def loci_from_breaks(self, ss_df, breaks_df):
        print('Determining loci from breaks')
        loci = pd.concat([ss_df.loc[(ss_df.chr == chr_n) & (ss_df.bp > start) & (ss_df.bp <= end), ['snp']].assign(locus=locus_n) for chr_n, chr_group in breaks_df.groupby('chr') for locus_n, (start, end) in chr_group[['start', 'end']].iterrows()])
        return(loci)
    
    def add_loci_from_breaks(self, sumstats, breaks):
        loci = self.loci_from_breaks(sumstats, breaks)
        return(sumstats.merge(loci, how='left', on='snp'))
    
    @staticmethod
    def consolidate(sets):
        '''
        function from: http://rosettacode.org/wiki/Set_consolidation#Python
        
        consolidates a list of sets. ex:
            [{1, 2}, {2}, {2, 3}, {2, 3}, {4}, {5, 6}, {6}]
            [{1, 2, 3}, {4}, {5, 6}]
            
        '''
        setlist = [s for s in sets if s]
        for i, s1 in enumerate(setlist):
            if s1:
                for s2 in setlist[i+1:]:
                    intersection = s1.intersection(s2)
                    if intersection:
                        s2.update(s1)
                        s1.clear()
                        s1 = s2
        return [s for s in setlist if s]
    
    def extend_loci(self, sumstats):
        overlaps = []
        for chr_n, chr_group in sumstats.groupby('chr'):
            for locus_n, locus_group in chr_group.groupby('locus'):
                lead_pos = locus_group.loc[locus_group.log10p.idxmax(), 'bp']
                
                in_window_idx = chr_group.index[(chr_group.bp - lead_pos).abs() <= self.delta]
                overlaps.append(set(sumstats.loc[in_window_idx, 'locus'].dropna().unique().astype(int)))
                sumstats.loc[in_window_idx, 'locus'] = locus_n
        
        sumstats['new_loci'] = np.nan            
        for i, overlap_loci in enumerate(self.consolidate(overlaps)):
            sumstats.loc[sumstats.locus.isin(overlap_loci), 'new_loci'] = i
        sumstats = sumstats.drop('locus', axis=1).rename({'new_loci':'locus'}, axis=1)
        return(sumstats)
    
    def prepare_gwas(self):
        # Load and clean
        sumstats = self.load_sumstats()
    
        # Get Loci
        breaks = self.load_breaks()
        sumstats = self.add_loci_from_breaks(sumstats, breaks)

        # Extend Loci
        sumstats = self.extend_loci(sumstats)
        sumstats = sumstats.dropna()
        
        # Save sumstats
        sumstats.to_csv(self.save_path, index=None, sep=self.ss_sep)
        return(sumstats)
    
#%%
if __name__ == '__main__':   
#%% Parser
    parser = argparse.ArgumentParser(description='Function to prepare chronic pain GWAS. If using with another GWAS, make sure assumptions hold')
    
    # Positional arguments
    parser.add_argument('--ss_path',
                        type=str,
                        help='path to raw summary statistics file')

    parser.add_argument('--breaks_path',
                        type=str,
                        help='path to loci breaks file') 
    
    parser.add_argument('--save_path',
                        type=str,
                        help='path to save final summary statistics file')
    
    # Optional arugments
    parser.add_argument('-ss_sep',
                        '--sumstats_separator',
                        type=str,
                        default='\t',
                        help='separator between summary statistics file columns')
    
    parser.add_argument('-breaks_sep',
                        '--breaks_separator',
                        type=str,
                        default='\t',
                        help='separator between loci breaks file columns')
    
    parser.add_argument('-delta',
                        '--window_size',
                        type=float,
                        default=0.5,
                        help='window size to extend break loci in Mb (5Kb -> 0.5Mb)')
    
    args = parser.parse_args()

    
    
#%% Main
    GWASPrep(ss_path=args.ss_path, breaks_path=args.breaks_path, save_path=args.save_path,
              ss_sep=args.sumstats_separator, breaks_sep=args.breaks_separator,
              delta=args.window_size/1e6).prepare_gwas()


