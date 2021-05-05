#%% libraries
from scipy.stats import norm
import numpy as np

# import matplotlib.pyplot as plt
#%% Function definitions
def p2z(p, betas):
    return(norm.isf(p/2)*np.sign(betas))
    
def z2p(z):
    return(norm.sf(np.abs(z))*2)

# def plot_loci(sumstats, ax=None):
#     if ax is None: ax = plt.gca()
#     sumstats.groupby(['chr', 'locus'])['snp'].nunique().unstack().plot(ax=ax, kind='bar', stacked=True, legend=False)