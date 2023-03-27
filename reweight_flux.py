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
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_7_reweight_flux'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_14_reweight_flux'
state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_reweight_flux'

weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_systs_truth_cuts'
#weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_3_systs_nuescat_cut'
#'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_30_systs_5bins'

sample_pot = 5.25536e17 #Size of sample used to make
nom_pot = 10e21 #Expected POT of sbnd
pot_sample_size = 1e9 #Renormalize x-axis to this value for nice numbers

pot_normalization = nom_pot/(sample_pot*pot_sample_size) #Get normalization to reach nominal flux

#Functions
def plot_reweighted_flux(universes,weights,**kwargs):
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
                r'Mean = ':f'{mean:.2f}',
                r'Std = ':f'{std:.2f}'}
  stats_unconstrained = plotters.convert_p_str(stats_unconstrained) #Get string version of stats
  stats_constrained = {f'After Constraint':'',
                r'Mean = ':f'{mean_reweighted:.2f}',
                r'Std = ':f'{std_reweighted:.2f}'}
  stats_constrained = plotters.convert_p_str(stats_constrained) #Get string version of stats

  fig,ax = plt.subplots(figsize=(6,6))
  ax.hist(universes,label='Unweighted',bins=bins,
          color='black',edgecolor='black',**kwargs)
  ax.hist(universes,weights=weights,label='Weighted',bins=bins,
          color='red',edgecolor='red',**kwargs)

  #Set labels
  trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
  ax.text(10,0.8,stats_unconstrained,transform=trans,fontsize=14,
          bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'))
  ax.text(10,0.62,stats_constrained,transform=trans,fontsize=14,
          bbox=dict(facecolor='none', edgecolor='red', boxstyle='round'))
  ax.set_xlabel(r'Flux($\nu/$m$^2/1\times10^9$ POT)')
  ax.set_ylabel('Prob. arb')
  ax.set_title(r'Constrainted Flux')
  plotters.set_style(ax)

  return fig,ax

def reweight_flux_enu(universes,edges,weights):
  """
  reweight flux for each energy bin by weight

  universes should be binned by energy

  Returns constrained and unconstrained bin heights, and fractional uncertainties
  """
  #Make arrays to fill
  height_constrained = np.zeros(len(edges)-1)
  height = height_constrained.copy()
  fracunc_constrained = height_constrained.copy()
  fracunc = height_constrained.copy()

  for i in range(height_constrained.shape[0]):
    flux_vals = universes[:,i] #Grab universe bin heights within energy bin
    height_constrained[i] = np.average(flux_vals, weights=weights)
    height[i] = np.average(flux_vals)
    fracunc_constrained[i] = np.sqrt(np.average((flux_vals - height_constrained[i]) ** 2, weights=weights))/height_constrained[i]
    fracunc[i] = np.sqrt(np.average((flux_vals - height[i]) ** 2))/height[i]

  return height_constrained,height,fracunc_constrained,fracunc

def plot_reweighted_flux_enu(height_constrained,height,edges,**kwargs):
  """
  Plot constrained and unconstrained flux values and their ratios
  """

  fig,(ax,ax_ratio) = plt.subplots(nrows=2,ncols=1,figsize=(6,6),gridspec_kw={'height_ratios': [2, 1]})
  ax.step(edges[:-1],height,color='black',label='Unconstrained',**kwargs)
  ax.step(edges[:-1],height_constrained,color='red',label='Constrained',**kwargs)
  ax.set_ylabel(r'Flux($\nu/$m$^2/1\times10^9$ POT/GeV)')
  ax.legend()

  ax_ratio.step(edges[:-1],height_constrained/height,color='black')
  ax_ratio.set_xlabel(r'True $E_\nu$ (GeV)')
  ax_ratio.set_ylabel(r'Con./Unc.')

  plotters.set_style(ax)
  plotters.set_style(ax_ratio)

  return fig,ax,ax_ratio


fname = 'state_all.root' #To get all universes and nominal one
print('Opening hist')
hist = uproot.open(f'{state_dir}/{fname}') #Open histograms
W_arr = np.loadtxt(f'{weights_dir}/weights.txt') #Get weights
W_arr_stat = np.loadtxt(f'{weights_dir}/weights_stat.txt') #Get weights
W_arr_syst = np.loadtxt(f'{weights_dir}/weights_syst.txt') #Get weights

#GEt bin values and edges
keys = hist.keys()
edges = hist[keys[0]].axis().edges() #They all have the same bin widths
values = np.zeros((len(keys),len(edges)-1)) #One less bin edge shape
for i,key in enumerate(keys): #Iterate over all universes
  values[i] = hist[key].values()

#Extract nominal and universe values
# nom = values[0] #I temporarily added the :-1 because the last bin is always zero
# #universes = np.concatenate((values[0:4],values[5:]),axis=0)
# universes = values[1:]
# n_universes = universes.shape[0]

nom = values[0] #I temporarily added the :-1 because the last bin is always zero
#universes = np.concatenate((values[0:4],values[5:]),axis=0)
universes = values[1:]
#edges = edges[:-1]
n_universes = universes.shape[0]

# plt.hist(np.sum(universes,axis=1)/(1e5*SBND_AREA)) #Flux normalized to area and 1e5 POT
# plt.savefig('tests/universe_flux.png')

save_dirs = ['total','stat_only','syst_only']
W_arrs = [W_arr,W_arr_stat,W_arr_syst]

for i,folder in enumerate(save_dirs):
  #print(f'{state_dir}/plots/{folder}')
  ws = W_arrs[i]
  #Get bin heights and uncertainties normalized to pot and area
  height_constrained,height,fracunc_constrained,fracunc = reweight_flux_enu(universes/(SBND_AREA)*pot_normalization,edges,ws)

  fig,ax = plot_reweighted_flux(np.sum(universes,axis=1)/(SBND_AREA)*pot_normalization,ws,histtype='step')
  plotters.save_plot('reweighted_flux',folder_name=f'{state_dir}/plots/{folder}',fig=fig)

  fig,ax,ax_ratio = plot_reweighted_flux_enu(height_constrained,height,edges)
  plotters.save_plot('reweighted_flux_enu',folder_name=f'{state_dir}/plots/{folder}',fig=fig)

  fig,ax,ax_ratio = plot_reweighted_flux_enu(fracunc_constrained,fracunc,edges)
  ax.set_ylabel('Fractional Uncertainty')
  plotters.save_plot('reweighted_fracunc_enu',folder_name=f'{state_dir}/plots/{folder}',fig=fig)

  plt.close('all')
