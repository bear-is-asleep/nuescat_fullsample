import pandas as pd
import seaborn as sns
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.transforms as transforms
from scipy.stats import chi2
import os

import sys
sys.path.append('../')
sys.path.append('/sbnd/app/users/brindenc/mypython/bc_utils')
from utils import plotters,pic
from CAFdata import *
import helpers
from datetime import date
import logging
logging.basicConfig(level=logging.CRITICAL)

day = date.today().strftime("%Y_%m_%d")
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_3_systs_nuescat_cut'
#'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_30_systs_5bins'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_systs_truth_cuts'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_29_systs_truth_cuts_eleeng'
state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_31_systs_truth_cuts_eleeng'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_4_3_systs_truth_cuts_eleeng_fv'
#Functions
def plot_reweighted_nue_dist(nom,universes,weights,**kwargs):
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

  fig,ax = plt.subplots(figsize=(8,6))
  ax.hist(universes,label='Unweighted',bins=bins,
          color='black',edgecolor='black',**kwargs)
  ax.hist(universes,weights=weights,label='Weighted',bins=bins,
          color='red',edgecolor='red',**kwargs)
  ax.axvline(nom,ls='--',color='blue')
  #ax.text(nom-50,60,f'CV = {nom:.0f}',fontsize=18)

  #Set labels
  trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
  ax.text(525,0.8,stats_unconstrained,transform=trans,fontsize=14,
          bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'))
  ax.text(525,0.62,stats_constrained,transform=trans,fontsize=14,
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
  xticks = ['']
  for i,edge in enumerate(edges):
    if i == len(edges)-1:break #Don't include last bin for out of bounds error
    xticks.append(f'{edge} - {edges[i+1]}')
  ax.set_xticklabels(xticks,rotation=30)
  ax.set_yticklabels(xticks,rotation=30)
  ax.set_xlabel(r'True $E_e$ [GeV]')
  ax.set_ylabel(r'True $E_e$ [GeV]')
  plt.colorbar()
  plotters.set_style(ax)
  return fig,ax

  
fname = 'state_Flux_multisim.root' #To get all universes and nominal one
hist = uproot.open(f'{state_dir}/{fname}') #Open histograms

#GEt bin values and edges
keys = hist.keys()
edges = hist[keys[0]].axis().edges() #They all have the same bin widths
values = np.zeros((len(keys),len(edges)-1)) #One less bin edge shape
for i,key in enumerate(keys): #Iterate over all universes
  values[i] = hist[key].values()

#Extract nominal and universe values
#print(np.sum(values,axis=1)[8])
data_ind = int(sys.argv[1]) #Universe number
meas = f'uni{data_ind}'
nom = values[0] #Nominal universe - all systs. are at mean value
data = values[data_ind] #Data we observe
universes = np.concatenate((values[0:data_ind],values[data_ind:]),axis=0)
#universes = values[1:]
#edges = edges[:-1]
n_universes = universes.shape[0]

plt.hist(np.sum(universes,axis=1))

#Calcu cov mat using Vij = 1/N sum_i((xu_i-xn_i)(xu_j-xn_j))
#covmat_stat = np.diag(np.mean(values,axis=0)) #stat
covmat_syst = np.zeros((len(nom),len(nom)))
covmat_stat = covmat_syst.copy()
mean_event_rates = np.mean(universes, axis=0)

for i,_ in enumerate(nom):
  for j,_ in enumerate(nom):     
    if i == j:
      covmat_stat[i,j] = np.mean((universes[:,i]))
    covmat_syst[i,j] = np.mean((universes[:,i]-mean_event_rates[i])*(universes[:,j]-mean_event_rates[j]))


#covmat = covmat_stat + covmat_syst #total covariance
#covmat = np.cov(values.T)
covmat = covmat_syst+covmat_stat
minerva_covmat = np.array([
  [98.7,1.22,1.72,1.38,0.42,-0.269],
  [1.22,27.3,1.63,1.14,.34,0.755],
  [1.72,1.63,40.1,1.88,0.596,1.35],
  [1.38,1.14,1.88,34.7,0.448,0.968],
  [0.42,0.34,0.596,0.448,18.9,0.778],
  [-0.269,0.755,1.35,0.968,0.778,59.5]])
minerva_edges = [0.8,2,3,5,7,9,20]
#covmat = np.identity(len(nom))
if True: #Save covariance matrices
  os.makedirs(f'{state_dir}/{meas}', exist_ok=True)
  np.savetxt(f'{state_dir}/{meas}/covmat_syst.txt',covmat_syst)
  np.savetxt(f'{state_dir}/{meas}/covmat_stat.txt',covmat_stat)
  np.savetxt(f'{state_dir}/{meas}/covmat.txt',covmat)
  np.savetxt(f'{state_dir}/minerva_covmat.txt',minerva_covmat)

  #Make plots
  fig,ax = plot_covmat(covmat_syst,edges,vmax=np.max(covmat))
  plotters.save_plot('covmat_syst',folder_name=f'{state_dir}/{meas}/plots/syst_only',fig=fig)

  fig,ax = plot_covmat(covmat_stat,edges,vmax=np.max(covmat))
  plotters.save_plot('covmat_stat',folder_name=f'{state_dir}/{meas}/plots/stat_only',fig=fig)

  fig,ax = plot_covmat(covmat,edges,vmax=np.max(covmat))
  plotters.save_plot('covmat',folder_name=f'{state_dir}/{meas}/plots/total',fig=fig)

  fig,ax = plot_covmat(minerva_covmat,minerva_edges,vmax=np.max(covmat))
  plotters.save_plot('minerva_covmat',folder_name=f'{state_dir}/{meas}/plots/total',fig=fig)
plt.close('all')
print('covmat made')

#Calc likelihood W = P(N_{\nu+e}|M) = \frac{1}{(2\pi)^{K/2}}\frac{1}{\sqrt{\Sigma_N}}
# \textrm{exp}\left( -\frac12 (\textbf{N}-\textbf{M})^T\Sigma_N^{-1} (\textbf{N}-\textbf{M})\right)
#χ2 =(N−M)T Σ−1(N−M)
def chi_weights_calc(covmat_input,data):
  dof = len(values[0]) #degrees of freedom
  chi_squared_arr = np.zeros(n_universes)
  normalization = 1/((2*np.pi)**(dof/2)*np.linalg.det(covmat_input)**(1/2))
  W_arr = np.zeros(n_universes)
  for i,uni in enumerate(universes):
    chi_squared_arr[i] = (data-uni).T@np.linalg.inv(covmat_input)@(data-uni)
    if i == 0: 
      #print('x2: ',chi_squared_arr[0])
      p_value = 1 - chi2.cdf(chi_squared_arr[0], dof)
    w = normalization*np.exp(-0.5*chi_squared_arr[i])
    W_arr[i] = w
  W_arr = len(universes)/np.sum(W_arr)*W_arr
  return W_arr,chi_squared_arr,p_value
W_arr,chi_squared_arr,p_value = chi_weights_calc(covmat,data)
#print('total: ',chi_squared_arr[0],p_value)
W_arr_syst,chi_squared_arr_syst,p_value_syst = chi_weights_calc(covmat_syst,data)
#print('syst: ',chi_squared_arr_syst[0], p_value_syst)
W_arr_stat,chi_squared_arr_stat,p_value_stat = chi_weights_calc(covmat_stat,data)
#print('stat: ',chi_squared_arr[0],p_value_stat)
file1 = open(f'{state_dir}/means.txt','a')
file1.writelines([str(d)+', ' for d in data])
file1.write('\n')
file1.close()

np.savetxt(f'{state_dir}/{meas}/data.txt',data)
np.savetxt(f'{state_dir}/{meas}/universes.txt',universes)
np.savetxt(f'{state_dir}/{meas}/x2.txt',chi_squared_arr)
np.savetxt(f'{state_dir}/{meas}/x2_stat.txt',chi_squared_arr_stat)
np.savetxt(f'{state_dir}/{meas}/x2_syst.txt',chi_squared_arr_syst)
#np.savetxt(f'{state_dir}/pvalue.txt',p_value)
#np.savetxt(f'{state_dir}/pvalue_stat.txt',p_value_stat)
#np.savetxt(f'{state_dir}/pvalue_syst.txt',p_value_syst)
np.savetxt(f'{state_dir}/{meas}/weights.txt',W_arr)
np.savetxt(f'{state_dir}/{meas}/weights_syst.txt',W_arr_syst)
np.savetxt(f'{state_dir}/{meas}/weights_stat.txt',W_arr_stat)
#print(W_arr,chi_squared_arr)
if True:
  plt.hist(W_arr)
  plotters.save_plot('ws',folder_name=f'{state_dir}/plots/total')
  plt.close()

  plt.hist(chi_squared_arr)
  plotters.save_plot('x2',folder_name=f'{state_dir}/plots/total')
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

  fig,ax = plot_reweighted_nue_dist(np.sum(data),
                                    np.sum(universes,axis=1),W_arr,histtype='step')
  plotters.save_plot('reweighted_ne',folder_name=f'{state_dir}/plots/total',fig=fig)

  fig,ax = plot_reweighted_nue_dist(np.sum(data),
                                    np.sum(universes,axis=1),W_arr_stat,histtype='step')
  plotters.save_plot('reweighted_ne',folder_name=f'{state_dir}/plots/stat_only',fig=fig)

  fig,ax = plot_reweighted_nue_dist(np.sum(data),
                                    np.sum(universes,axis=1),W_arr_syst,histtype='step')
  plotters.save_plot('reweighted_ne',folder_name=f'{state_dir}/plots/syst_only',fig=fig)