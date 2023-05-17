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
import logging
logging.basicConfig(level=logging.CRITICAL)

SBND_AREA = 1

day = date.today().strftime("%Y_%m_%d")
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_7_reweight_flux'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_14_reweight_flux'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_reweight_flux'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_29_reweight_flux_stride100'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_30_reweight_flux_NCE'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_30_reweight_flux_var'
#state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_31_reweight_flux_var'
state_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_31_reweight_all_flux_var'

#weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_systs_truth_cuts'
#weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_3_systs_nuescat_cut'
#weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_29_systs_truth_cuts_eleeng'
weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_31_systs_truth_cuts_eleeng'
#weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_4_3_systs_truth_cuts_eleeng_fv'
#'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_30_systs_5bins'

#sample_pot = 1.55131e+19 #Size of sample used to make
#sample_pot = 6.16911e+19
sample_pot = 3.08e20 #Total from MCP2022A
nom_pot = 10e21 #Expected POT of sbnd
pot_sample_size = 1e9 #Renormalize x-axis to this value for nice numbers

#ratio_conv = 394766.0450858611 #This comes from my conversion from gev to cm being incorrect when generating the sample
mm2_to_cm2 = 1e2 #I converted to mm instead of cm
pot_normalization = mm2_to_cm2/sample_pot #Get normalization to reach nominal flux

#Functions
def plot_reweighted_flux(nom,universes,weights,**kwargs):
  """
  plot number of nu e scatters for all universes and weights

  make sure to sum across all bins so universes is (uni) not (uni,bin) dim
  """
  plot_normalization = 1e9
  universes*=plot_normalization
  nom*=plot_normalization

  mean = np.average(universes)
  std = np.std(universes)

  mean_reweighted = np.average(universes, weights=weights)
  std_reweighted = np.sqrt(np.average((universes - mean_reweighted) ** 2, weights=weights))

  bw = std/3 #bin width
  bins = np.arange(np.min(universes-bw),np.max(universes)+2*bw,bw)

  stats_unconstrained = {f'Before Constraint':'',
                r'Mean = ':f'{mean:.2f}',
                r'Std = ':f'{std:.2f}'}
  stats_unconstrained = plotters.convert_p_str(stats_unconstrained) #Get string version of stats
  stats_constrained = {f'After Constraint':'',
                r'Mean = ':f'{mean_reweighted:.2f}',
                r'Std = ':f'{std_reweighted:.2f}'}
  stats_constrained = plotters.convert_p_str(stats_constrained) #Get string version of stats

  #Write mean  value to txt file
  file1 = open(f'{state_dir}/means.txt','a')
  file1.write(f'{mean_reweighted}, {std_reweighted} \n')
  file1.close()

  fig,ax = plt.subplots(figsize=(6,6))
  ax.hist(universes,label='Unweighted',bins=bins,
          color='black',edgecolor='black',**kwargs)
  ax.hist(universes,weights=weights,label='Weighted',bins=bins,
          color='red',edgecolor='red',**kwargs)
  ax.axvline(nom,ls='--',label='Data')
  #Set labels
  trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
  ax.text(6.6e0,0.8,stats_unconstrained,transform=trans,fontsize=14,
          bbox=dict(facecolor='none', edgecolor='black', boxstyle='round'))
  ax.text(6.6e0,0.62,stats_constrained,transform=trans,fontsize=14,
          bbox=dict(facecolor='none', edgecolor='red', boxstyle='round'))
  ax.set_xlabel(r'Flux($\nu_\mu/$cm$^2$/ 10$^{9}$POT)')
  ax.set_ylabel('Prob. arb')
  ax.set_title(r'Constrainted Flux')
  plotters.set_style(ax)

  return fig,ax

def reweight_flux_enu(nom,universes,edges,weights):
  """
  reweight flux for each energy bin by weight

  universes should be binned by energy

  Returns constrained and unconstrained bin heights, and fractional uncertainties
  """
  #Make arrays to fill
  height_constrained = np.zeros(len(edges)-1)
  height = height_constrained.copy()

  #Resolution
  fracunc_constrained = height_constrained.copy()
  fracunc = height_constrained.copy()

  #Bias
  fracbias_constrained = height_constrained.copy()
  fracbias = height_constrained.copy()

  for i in range(height_constrained.shape[0]):
    flux_vals = universes[:,i] #Grab universe bin heights within energy bin
    flux_nom = nom[i] #Central value flux
    #Get heights
    height_constrained[i] = np.average(flux_vals, weights=weights)
    height[i] = np.average(flux_vals)
    #Get resolutions
    fracunc_constrained[i] = np.sqrt(np.average((flux_vals - height_constrained[i]) ** 2, weights=weights))/height_constrained[i]
    fracunc[i] = np.sqrt(np.average((flux_vals - height[i]) ** 2))/height[i]
    #Get biases
    fracbias_constrained[i] = np.sqrt(np.average((flux_vals - flux_nom) ** 2,weights=weights))/flux_nom
    fracbias[i] = np.sqrt(np.average((flux_vals - flux_nom) ** 2))/flux_nom

  return height_constrained,height,fracunc_constrained,fracunc,fracbias_constrained,fracbias

def plot_reweighted_flux_enu(nom,height_constrained,height,edges,fracunc_constrained=None,fracunc=None,**kwargs):
  """
  Plot constrained and unconstrained flux values and their ratios
  """
  # Calculate the central values of each bin
  bin_widths = np.diff(edges)
  central_values = edges[:-1] - bin_widths/2

  fig,(ax,ax_ratio) = plt.subplots(nrows=2,ncols=1,figsize=(10,8),gridspec_kw={'height_ratios': [2, 1]})
  ax.step(edges[:-1],height,color='black',label='Unconstrained',**kwargs)
  ax.step(edges[:-1],height_constrained,color='red',label='Constrained',**kwargs)
  if nom is not None:
    ax.step(edges[:-1],nom,ls='--',color='blue',label='Data')
  if fracunc_constrained is not None:
    y_err = fracunc_constrained*height_constrained 
    ax.errorbar(central_values, height_constrained, yerr=y_err, capsize=2,fmt='none',color='red',alpha=0.3)
  if fracunc is not None:
    y_err = fracunc*height 
    ax.errorbar(central_values, height, yerr=y_err, capsize=2,fmt='none',color='black',alpha=0.3)
    
  ax.set_ylabel(r'Flux($\nu_\mu/$cm$^2$/POT)')
  ax.legend()

  ax_ratio.step(edges[:-1],height_constrained/height,color='black')
  ax_ratio.set_xlabel(r'True $E_\nu$ (GeV)')
  ax_ratio.set_ylabel(r'Con./Unc.')

  plotters.set_style(ax)
  plotters.set_style(ax_ratio)

  return fig,ax,ax_ratio

def plot_reweighted_flux_resolution_bias(nom,height_constrained,height,edges,fracunc_constrained=None,fracunc=None,
                                         fracbias_constrained=None,fracbias=None,**kwargs):
  """
  Plot constrained and unconstrained flux values and their ratios
  """
  # Calculate the central values of each bin
  bin_widths = np.diff(edges)
  central_values = edges[:-1] - bin_widths/2

  fig,(ax,ax_unc_bias) = plt.subplots(nrows=2,ncols=1,figsize=(10,8),gridspec_kw={'height_ratios': [1, 1]})
  ax.step(edges[:-1],height,color='black',label='Unconstrained',**kwargs)
  if height_constrained is not None:
    ax.step(edges[:-1],height_constrained,color='red',label='Constrained',**kwargs)
  if nom is not None:
    ax.step(edges[:-1],nom,color='blue',ls='--',label='Data')
  if fracunc_constrained is not None and fracbias_constrained is not None:
    y_err = np.sqrt(fracunc_constrained**2+fracbias_constrained**2)*height_constrained 
    ax.errorbar(central_values, height_constrained, yerr=y_err, capsize=2,fmt='none',color='red',alpha=0.3)
  if fracunc is not None and fracbias is not None:
    y_err = np.sqrt(fracunc**2+fracbias**2)*height 
    ax.errorbar(central_values, height, yerr=y_err, capsize=2,fmt='none',color='black',alpha=0.3)
    
  ax.set_ylabel(r'Flux($\nu_\mu/$cm$^2$/POT)')
  ax.legend()

  if fracunc is not None:
    ax_unc_bias.step(edges[:-1],fracunc,color='black',ls='--',label='Resolution')
  if fracbias is not None:
    ax_unc_bias.step(edges[:-1],fracbias,color='black',ls='-',label='Bias')
  if fracunc_constrained is not None:
    ax_unc_bias.step(edges[:-1],fracunc_constrained,color='red',ls='--',label='Constrained Resolution')
  if fracbias_constrained is not None:
    ax_unc_bias.step(edges[:-1],fracbias_constrained,color='red',ls='-',label='Constrained Bias')
  ax_unc_bias.set_xlabel(r'True $E_\nu$ (GeV)')
  ax_unc_bias.set_ylabel(r'Fractional Uncertainty')
  ax_unc_bias.legend()

  plotters.set_style(ax)
  plotters.set_style(ax_unc_bias,legend_size=12)

  return fig,ax,ax_unc_bias

def plot_reweighted_flux_total_unc(nom,height_constrained,height,edges,fracunc_constrained=None,fracunc=None,
                                         fracbias_constrained=None,fracbias=None,**kwargs):
  """
  Plot constrained and unconstrained flux values and their ratios
  """
  # Calculate the central values of each bin
  bin_widths = np.diff(edges)
  central_values = edges[:-1] - bin_widths/2

  fig,(ax,ax_unc_bias) = plt.subplots(nrows=2,ncols=1,figsize=(10,8),gridspec_kw={'height_ratios': [2, 1]})
  ax.step(edges[:-1],height,color='black',label='Unconstrained',**kwargs)
  if height_constrained is not None:
    ax.step(edges[:-1],height_constrained,color='red',label='Constrained',**kwargs)
  if nom is not None:
    ax.step(edges[:-1],nom,color='blue',ls='--',label='Data')
  if fracunc_constrained is not None and fracbias_constrained is not None:
    y_err = np.sqrt(fracunc_constrained**2+fracbias_constrained**2)*height_constrained 
    ax.errorbar(central_values, height_constrained, yerr=y_err, capsize=2,fmt='none',color='red',alpha=0.3)
  if fracunc is not None and fracbias is not None:
    y_err = np.sqrt(fracunc**2+fracbias**2)*height 
    ax.errorbar(central_values, height, yerr=y_err, capsize=2,fmt='none',color='black',alpha=0.3)
    
  ax.set_ylabel(r'Flux($\nu_\mu/$cm$^2$/POT)')
  ax.legend()

  if fracunc is not None and fracbias is not None:
    ax_unc_bias.step(edges[:-1],np.sqrt(fracunc**2+fracbias**2),color='black',label='Unconstrained')
  if fracunc_constrained is not None and fracbias_constrained is not None:
    ax_unc_bias.step(edges[:-1],np.sqrt(fracunc_constrained**2+fracbias_constrained**2),color='red',label='Constrained')
  ax_unc_bias.set_xlabel(r'True $E_\nu$ (GeV)')
  ax_unc_bias.set_ylabel(r'Fractional Uncertainty')
  ax_unc_bias.legend()

  plotters.set_style(ax)
  plotters.set_style(ax_unc_bias,legend_size=12)

  return fig,ax,ax_unc_bias

data_ind = int(sys.argv[1]) #data ind

fname = 'state_total_flux.root' #To get all universes and nominal one
print('Opening hist')
hist = uproot.open(f'{state_dir}/{fname}') #Open histograms
#meas = 'uni8_numpy_covmat/' #What data to use
meas = f'uni{data_ind}'
W_arr = np.loadtxt(f'{weights_dir}/{meas}/weights.txt') #Get weights
W_arr_stat = np.loadtxt(f'{weights_dir}/{meas}/weights_stat.txt') #Get weights
W_arr_syst = np.loadtxt(f'{weights_dir}/{meas}/weights_syst.txt') #Get weights

#GEt bin values and edges
keys = hist.keys()
edges = hist[keys[0]].axis().edges() #They all have the same bin widths
values = np.zeros((len(keys),len(edges)-1)) #One less bin edge shape
for i,key in enumerate(keys): #Iterate over all universes
  #print(hist[key].values()[1],',')
  values[i] = hist[key].values()

#Extract nominal and universe values
data = values[data_ind] #I temporarily added the :-1 because the last bin is always zero
universes = np.concatenate((values[0:data_ind],values[data_ind:]),axis=0)
universes = values[1:]
nom = values[0]
#nom = data
# n_universes = universes.shape[0]
#Iterations*2 bin counts and modify edges
iterations = 0
if iterations != 0:
  new_values = np.zeros((values.shape[0],(int(values.shape[1]/(2**iterations))+values.shape[1]%2)))
  for i,_ in enumerate(values):
    new_values[i],new_edges = helpers.combine_bins_sum(values[i],edges,iterations=iterations)
  values = new_values
  edges=new_edges
#nom = values[0] #I temporarily added the :-1 because the last bin is always zero
#universes = np.concatenate((values[0:4],values[5:]),axis=0)
#universes = values[1:]
#edges = edges[:-1]
n_universes = universes.shape[0]

# plt.hist(np.sum(universes,axis=1)/(1e5*SBND_AREA)) #Flux normalized to area and 1e5 POT
# plt.savefig('tests/universe_flux.png')

#save_dirs = ['total','stat_only','syst_only']
save_dirs = ['stat_only']
W_arrs = [W_arr[:100],W_arr_stat[:100],W_arr_syst[:100]]

nom = nom*pot_normalization
data = data*pot_normalization
universes = universes*pot_normalization

for i,folder in enumerate(save_dirs):
  #print(f'{state_dir}/plots/{folder}')
  ws = W_arrs[i]
  #Get bin heights and uncertainties normalized to pot and area
  height_constrained,height,fracunc_constrained,fracunc,fracbias_constrained,fracbias= reweight_flux_enu(nom,universes,edges,ws)
  fig,ax = plot_reweighted_flux(np.sum(data),
                                np.sum(universes,axis=1),ws,histtype='step')
  plotters.save_plot(f'reweighted_flux_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)

  fig,ax,ax_ratio = plot_reweighted_flux_enu(data,height_constrained,height,edges,
                                             fracunc_constrained=fracunc_constrained,fracunc=fracunc)
  ax.set_xlim([0,3])
  ax_ratio.set_xlim([0,3])
  plotters.save_plot(f'reweighted_flux_enu_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)
  ax.set_yscale('log')
  plotters.save_plot(f'reweighted_flux_enu_logy_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)

  fig,ax,ax_ratio = plot_reweighted_flux_enu(None,np.sqrt(fracunc_constrained**2+fracbias_constrained**2),
                                             np.sqrt(fracunc**2+fracbias**2),edges)
  ax.set_ylabel('Fractional Uncertainty')
  ax.set_xlim([0,3])
  ax_ratio.set_xlim([0,3])
  plotters.save_plot(f'reweighted_fracunc_enu_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)

  fig,ax,ax_unc = plot_reweighted_flux_resolution_bias(data,height_constrained,height,edges,fracunc_constrained,fracunc,
                                         fracbias_constrained,fracbias)
  ax.set_xlim([0,3])
  ax_unc.set_xlim([0,3])
  plotters.save_plot(f'reweighted_flux_fracunc_constrained_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)
  ax.set_yscale('log')
  plotters.save_plot(f'reweighted_flux_fracunc_constrained_logy_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)

  fig,ax,ax_unc = plot_reweighted_flux_resolution_bias(data,height_constrained,height,edges,fracunc_constrained,fracunc,
                                         fracbias_constrained,fracbias)
  ax.set_xlim([0,3])
  ax_unc.set_xlim([0,3])
  plotters.save_plot(f'reweighted_flux_fracunc_constrained_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)
  ax.set_yscale('log')
  plotters.save_plot(f'reweighted_flux_fracunc_constrained_logy_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)

  fig,ax,ax_unc = plot_reweighted_flux_total_unc(data,height_constrained,height,edges,fracunc_constrained,fracunc,
                                         fracbias_constrained,fracbias)
  ax.set_xlim([0,3])
  ax_unc.set_xlim([0,3])
  plotters.save_plot(f'reweighted_flux_fracunc_total_constrained_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)
  ax.set_yscale('log')
  plotters.save_plot(f'reweighted_flux_fracunc_total_constrained_logy_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)
  
  #No bias
  fig,ax,ax_ratio = plot_reweighted_flux_enu(None,fracunc_constrained,
                                             fracunc,edges)
  ax.set_title('Uncertainty w/ no Bias')
  ax.set_ylabel('Fractional Uncertainty')
  ax.set_xlim([0,3])
  ax_ratio.set_xlim([0,3])
  plotters.save_plot(f'reweighted_fracunc_enu_nobias_{meas}',folder_name=f'{state_dir}/plots/{folder}',fig=fig)
  
  
  plt.close('all')
