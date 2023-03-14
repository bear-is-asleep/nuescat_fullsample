import pandas as pd
import seaborn as sns
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.transforms as transforms

import sys
sys.path.append('../')
sys.path.append('/sbnd/app/users/brindenc/mypython/bc_utils')
from utils import plotters,pic
from CAFdata import *
import helpers
from datetime import date

day = date.today().strftime("%Y_%m_%d")
state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_3_systs_nuescat_cut'
#'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_30_systs_5bins'

#Functions
def plot_reweighted_nue_dist(universes,weights,**kwargs):
  """
  plot number of nu e scatters for all universes and weights

  make sure to sum across all bins so universes is (uni) not (uni,bin) dim
  """
  mean = np.average(universes)
  std = np.std(universes)

  mean_reweighted = np.average(universes, weights=weights)
  std_reweighted = np.sqrt(np.average((universes - mean_reweighted) ** 2, weights=weights))

  bw = std/5 #bin width
  bins = np.arange(np.min(universes-bw),np.max(universes)+2*bw,bw)

  stats_unconstrained = {f'Before Constraint':'',
                r'$N_{\nu+e}$ Mean = ':f'{mean:.1f}',
                r'$N_{\nu+e}$ Std = ':f'{std:.1f}'}
  stats_unconstrained = plotters.convert_p_str(stats_unconstrained) #Get string version of stats
  stats_constrained = {f'After Constraint':'',
                r'$N_{\nu+e}$ Mean = ':f'{mean_reweighted:.1f}',
                r'$N_{\nu+e}$ Std = ':f'{std_reweighted:.1f}'}
  stats_constrained = plotters.convert_p_str(stats_constrained) #Get string version of stats

  fig,ax = plt.subplots(figsize=(6,6))
  ax.hist(universes,label='Unweighted',bins=bins,
          color='black',edgecolor='black',**kwargs)
  ax.hist(universes,weights=W_arr,label='Weighted',bins=bins,
          color='red',edgecolor='red',**kwargs)

  #Set labels
  trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
  ax.text(90,0.8,stats_unconstrained,transform=trans,fontsize=14,
          bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'))
  ax.text(90,0.62,stats_constrained,transform=trans,fontsize=14,
          bbox=dict(facecolor='none', edgecolor='red', boxstyle='round'))
  ax.set_xlabel(r'$N_{\nu+e}$')
  ax.set_ylabel('Prob. arb')
  ax.set_title(r'$N_{\nu+e}$ Constraint')
  plotters.set_style(ax)

  return fig,ax

def plot_covmat(covmat,edges,**kwargs):
  """
  Plot covmat
  """
  fig,ax = plt.subplots(figsize=(8,8),tight_layout=True)
  plt.imshow(covmat,cmap='viridis',**kwargs)
  # add text to display the value of each bin
  for i in range(covmat.shape[0]):
    for j in range(covmat.shape[1]):
      plt.text(j, i, f'{covmat[i,j]:.1f}', ha='center', va='center', color='white')
  plt.xticks(np.arange(len(edges)-1),edges[:-1]) #Label ticks with energy spectrum
  plt.yticks(np.arange(len(edges)-1),edges[:-1]) #Label ticks with energy spectrum
  plt.colorbar()
  return fig,ax


  
fname = 'state_all.root' #To get all universes and nominal one
hist = uproot.open(f'{state_dir}/{fname}') #Open histograms

#GEt bin values and edges
keys = hist.keys()
edges = hist[keys[0]].axis().edges() #They all have the same bin widths
values = np.zeros((len(keys),len(edges)-1)) #One less bin edge shape
for i,key in enumerate(keys): #Iterate over all universes
  values[i] = hist[key].values()

#Extract nominal and universe values
nom = values[0,:-1] #I temporarily added the :-1 because the last bin is always zero
universes = values[1:,:-1]
edges = edges[:-1]
n_universes = universes.shape[0]

plt.hist(np.sum(universes,axis=1))
plt.savefig('tests/universe.png')

#Calcu cov mat using Vij = 1/N sum_i((xu_i-xn_i)(xu_j-xn_j))
covmat_syst = np.zeros((len(nom),len(nom))) #syst
covmat_stat = np.zeros((len(nom),len(nom))) #stat

for _,uni in enumerate(universes):
  for i,_ in enumerate(nom):
    for j,_ in enumerate(nom):
      if (uni == nom).all(): continue #Don't include nominal universe in calc
      covmat_syst[i][j] += (uni[i]-nom[i])*(uni[j]-nom[j])
      if i == j:
        covmat_stat[i][j] += np.sqrt(uni[i])

covmat_syst/=n_universes
covmat_stat/=n_universes
covmat = covmat_stat+covmat_syst #total covariance
if True: #Save covariance matrices
  np.savetxt(f'{state_dir}/covmat_syst.txt',covmat_syst)
  np.savetxt(f'{state_dir}/covmat_stat.txt',covmat_stat)
  np.savetxt(f'{state_dir}/covmat.txt',covmat)

  #Make plots
  fig,ax = plot_covmat(covmat_syst,edges,vmax=np.max(covmat))
  plotters.save_plot('covmat_syst',folder_name=f'{state_dir}/plots',fig=fig)

  fig,ax = plot_covmat(covmat_stat,edges,vmax=np.max(covmat))
  plotters.save_plot('covmat_stat',folder_name=f'{state_dir}/plots',fig=fig)

  fig,ax = plot_covmat(covmat,edges,vmax=np.max(covmat))
  plotters.save_plot('covmat',folder_name=f'{state_dir}/plots',fig=fig)

print('covmat made')

#Calc likelihood W = P(N_{\nu+e}|M) = \frac{1}{(2\pi)^{K/2}}\frac{1}{\sqrt{\Sigma_N}}
# \textrm{exp}\left( -\frac12 (\textbf{N}-\textbf{M})^T\Sigma_N^{-1} (\textbf{N}-\textbf{M})\right)
#χ2 =(N−M)T Σ−1(N−M)
dof = len(values[0]) #degrees of freedom
chi_squared_arr = np.zeros(n_universes)
normalization = 1/((2*np.pi)**(dof/2)*np.linalg.det(covmat)**(1/2))
W_arr = np.zeros(n_universes)
#universes_reweighted = np.zeros(universes.shape)
norm_bayes = np.zeros(universes.shape[1]) #int pi(m')p(n|m')dm' = sum_i u_i/nom (vector) p(n_i|m)
for i,uni in enumerate(universes):
  chi_squared_arr[i] = (nom-uni).T@np.linalg.inv(covmat)@(nom-uni)
  w = normalization*np.exp(-0.5*chi_squared_arr[i])
  W_arr[i] = w
  #universes_reweighted[i] = w*uni
  #norm_bayes += universes_reweighted[i]/nom
plt.close('all')
#for i,uni in enumerate(universes_reweighted):
#  print(uni/norm_bayes,uni,norm_bayes)
#universes_reweighted = [uni/norm_bayes for uni in universes_reweighted] #normalize
#universes_reweighted = np.sum(universes,axis=0)/np.sum(universes_reweighted,axis=0)*universes_reweighted #normalize
#np.savetxt('tests/reweighted.txt',universes_reweighted)
W_arr = len(universes)/np.sum(W_arr)*W_arr
np.savetxt(f'{state_dir}/universes.txt',universes)
np.savetxt(f'{state_dir}/x2.txt',np.exp(-0.5*chi_squared_arr))
np.savetxt(f'{state_dir}/weights.txt',W_arr)
#print(W_arr,chi_squared_arr)
if True:
  plt.hist(W_arr)
  plotters.save_plot('ws',folder_name=f'{state_dir}/plots')
  plt.close()

  plt.hist(chi_squared_arr)
  plotters.save_plot('x2',folder_name=f'{state_dir}/plots')
  plt.close()

  
  plt.hist(np.sum(universes,axis=1),
           alpha=0.5,
           label='Unweighted')
  plt.hist(np.sum(universes,axis=1),
           weights=W_arr,
           alpha=0.5,
           label='Weighted-nonnormalized')
  plt.legend()
  plotters.save_plot('n_ele',folder_name=f'{state_dir}/plots')
  plt.close()

  fig,ax = plot_reweighted_nue_dist(np.sum(universes,axis=1),W_arr,histtype='step')
  plotters.save_plot('reweighted_ne',folder_name=f'{state_dir}/plots',fig=fig)