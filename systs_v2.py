import numpy as np
import uproot
import pandas as pd
import matplotlib.pyplot as plt
import helpers
import os
import matplotlib.cm as cm
plt.style.use(['science','no-latex'])

import sys
sys.path.append('/sbnd/app/users/brindenc/mypython/bc_utils')
from utils import plotters,pic

sample = ''

day = helpers.day
sample_pot = 3.08e20
mm2_to_cm2 = 1e2 #I converted to mm instead of cm
pot_normalization = mm2_to_cm2/sample_pot #Get normalization to reach nominal flux
#suffix = '_systs_nuescat_cut'
#state_folder = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_3_systs_nuescat_cut'
#state_folder = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_systs_truth_cuts'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_reweight_flux'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_29_systs_truth_cuts_eleeng'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_31_systs_truth_cuts_eleeng'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_29_reweight_flux_stride100'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_31_reweight_all_flux_var'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_4_5_systs_truth_cuts_nueng_av'
state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_4_6_systs_truth_cuts_nueng_av'

weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_31_systs_truth_cuts_eleeng/cv'
W_arr = np.loadtxt(f'{weights_dir}/weights.txt') #Get weights
W_arr_stat = np.loadtxt(f'{weights_dir}/weights_stat.txt') #Get weights
W_arr_syst = np.loadtxt(f'{weights_dir}/weights_syst.txt') #Get weights

W_arr_stat = W_arr_stat[:1000]

plot_folder = f'{state_folder}/plots/systs/'

def sqrt_sum_of_squares(numbers):
  """
  Takes a list of numbers and returns the sum of their squares.
  """
  return np.sqrt(np.sum([x**2 for x in numbers]))

def plot_resolution_bias_stat(nom_heights,heights,edges,labels,colors,heights_constrained=None,fracunc=None,fracbias=None,fracunc_constrained=None,
                                         fracbias_constrained=None,fracstat=None,fracstat_constrained=None,**kwargs):
  """
  Plot constrained and unconstrained values and their ratios
  """

  # Calculate the central values of each bin
  bin_widths = np.diff(edges)
  central_values = edges[1:] - bin_widths/2

  fig,(ax,ax_unc_bias) = plt.subplots(nrows=2,ncols=1,figsize=(10,8),gridspec_kw={'height_ratios': [1, 1]})
  for i,_ in enumerate(heights):
    # ax.step(central_values,nom_heights,where='mid',
    #         color=colors(i),label=labels[i],**kwargs)
    if fracunc_constrained is not None:
      y_err = np.sqrt(fracunc_constrained**2+fracbias_constrained**2)*heights_constrained 
      ax.errorbar(central_values, heights_constrained, yerr=y_err, capsize=2,fmt='none',color='red',alpha=0.3)
    if fracunc is not None:
      if fracbias is not None:
        #print(fracbias.shape,fracunc.shape)
        fracunc[i] = np.sqrt(fracbias[i]**2+fracunc[i]**2)
      #total_unc = np.nansum(fracunc[i]*heights[-1]*100/np.sum(heights[-1]))
      total_unc = np.nansum(fracunc[i])*100
      #print(total_unc,fracunc[i],heights[-1])
      label = f'{labels[i]} ({total_unc:.1f}%)'
      y_err = fracunc[i]*nom_heights 
      #ax.errorbar(central_values, nom_heights,yerr=y_err,capsize=2,fmt='none',color=colors(i),alpha=0.5)
      if labels[i] == 'Total':
        ax_unc_bias.step(edges[:-1],fracuncs[i],where='post',color='black',linewidth=2,ls='-',label=label)
        ax_unc_bias.step([edges[-2], edges[-1]], [fracuncs[i,-1],fracuncs[i,-1]],where='post',color='black',linewidth=2,ls='-')
      else:
        ax_unc_bias.step(edges[:-1],fracuncs[i],where='post',color=colors(i),linewidth=1,ls='-',label=label)
        ax_unc_bias.step([edges[-2], edges[-1]], [fracuncs[i,-1],fracuncs[i,-1]],where='post',color=colors(i),linewidth=1,ls='-')
  if nom_heights is not None:
    y_err = fracunc[i]*heights[-1] 
    ax.step(edges[:-1],nom_heights,where='post',color='blue',label='Nominal',
            linewidth=3,**kwargs)
    #Add the last line
    ax.step([edges[-2], edges[-1]], [nom_heights[-1],nom_heights[-1]], where='post',color='blue',
            linewidth=3,**kwargs)
    ax.errorbar(central_values, heights[-1],yerr=y_err,capsize=2,fmt='none',color='black')

  ax.set_ylabel(r'Flux($\nu_\mu/$cm$^2$/POT)')
  #ax.set_ylabel(r'$N_{\nu+e}$')
  ax.legend()

  ax_unc_bias.set_xlabel(r'True $E_e$ (GeV)')
  ax_unc_bias.set_ylabel(r'Fractional Uncertainty')
  ax_unc_bias.legend()

  plotters.set_style(ax)
  plotters.set_style(ax_unc_bias,legend_size=12)

  return fig,ax,ax_unc_bias

def combine_frac_height_vals(indeces,heights,heights_constrained,fracuncs_constrained,fracuncs,fracbiases_constrained,
                      fracbiases,fracstats_constrained,fracstats,return_combined=True):
  """
  Combine fractional uncertainty values and return the combined values at the end of the list
  """
  #U: universes
  #E: edges - 1

  #Uncertainties to combine
  heights_constrained_combined = np.zeros(heights[0].shape) #E
  heights_combined = heights_constrained_combined.copy()
  fracuncs_constrained_combined = np.zeros(fracbiases[0].shape) #E
  fracuncs_combined = fracuncs_constrained_combined.copy()
  fracbiases_constrained_combined = fracuncs_constrained_combined.copy()
  fracbiases_combined = fracuncs_constrained_combined.copy()
  fracstats_constrained_combined = fracuncs_constrained_combined.copy()
  fracstats_combined = fracuncs_constrained_combined.copy()

  #Combine
  for iuni,_ in enumerate(fracuncs_constrained_combined):
    for iedge,_ in enumerate(fracuncs_constrained_combined):
      #Compute average between heights
      heights_constrained_combined[iedge] = np.average([heights_constrained[ind][iedge] for ind in indeces])
      heights_combined[iedge] = np.average([heights[ind][iedge] for ind in indeces])
      #Compute the combined fractional uncertainty
      fracuncs_constrained_combined[iedge] = sqrt_sum_of_squares([fracuncs_constrained[ind][iedge] for ind in indeces])
      fracuncs_combined[iedge] = sqrt_sum_of_squares([fracuncs[ind][iedge] for ind in indeces])
      fracbiases_constrained_combined[iedge] = sqrt_sum_of_squares([fracbiases_constrained[ind][iedge] for ind in indeces])
      fracbiases_combined[iedge] = sqrt_sum_of_squares([fracbiases[ind][iedge] for ind in indeces])
      fracstats_constrained_combined[iedge] = sqrt_sum_of_squares([fracstats_constrained[ind][iedge] for ind in indeces])
      fracstats_combined[iedge] = sqrt_sum_of_squares([fracstats[ind][iedge] for ind in indeces])
  if return_combined:
    return heights_constrained_combined,heights_combined,fracuncs_constrained_combined,fracuncs_combined,fracbiases_constrained_combined,fracbiases_combined,fracstats_constrained_combined,fracstats_combined
  else:
    #Drop the uncertainties we just combined
    heights_constrained = helpers.drop_matching_indices(heights_constrained,indeces)
    heights = helpers.drop_matching_indices(heights,indeces)
    fracuncs_constrained = helpers.drop_matching_indices(fracuncs_constrained,indeces)
    fracuncs = helpers.drop_matching_indices(fracuncs,indeces)
    fracbiases_constrained = helpers.drop_matching_indices(fracbiases_constrained,indeces)
    fracbiases = helpers.drop_matching_indices(fracbiases,indeces)
    fracstats_constrained = helpers.drop_matching_indices(fracstats_constrained,indeces)
    fracstats = helpers.drop_matching_indices(fracstats,indeces)

    #Append combined uncertainties to the end of the list
    heights_constrained.append(heights_constrained_combined)
    heights.append(heights_combined)
    fracuncs_constrained.append(fracuncs_constrained_combined)
    fracuncs.append(fracuncs_combined)
    fracbiases_constrained.append(fracbiases_constrained_combined)
    fracbiases.append(fracbiases_combined)
    fracstats_constrained.append(fracstats_constrained_combined)
    fracstats.append(fracstats_combined)
    return heights_constrained,heights,fracuncs_constrained,fracuncs,fracbiases_constrained,fracbiases,fracstats_constrained,fracstats

fnames = os.listdir(state_folder)
state_names = [] #Label names
state_fnames = [] #File names
for i,fname in enumerate(fnames):
  if fname == 'state_all.root': continue #skip
  if fname[:5] == 'state' and fname[-5:] == '.root' and 'band' not in fname and 'nom' not in fname: #Get state files
    state_names.append(fname[6:-5])
    state_fnames.append(fname)
print(state_fnames)

#Store information for plotting
heights_constrained = []
heights = []
fracuncs_constrained = []
fracuncs = []
fracbiases_constrained = []
fracbiases = []
fracstats_constrained = []
fracstats = []

#Indeces to drop,keep total indeces, etc.
optics_indeces = [0,1]
kplus_index = [3]
piplus_index = [12]
kflux_indeces = [2,4]
xsec_indeces = [7,11]
piminus_indeces = [8]
other_indeces = xsec_indeces + piminus_indeces + kflux_indeces
flux_index = [15]
drop_indeces = [5,6,9,10,-1] #Components of XSec and unweighted

og_labels = ['expskin_Flux', 
          'horncurrent_Flux', 
          'kminus_Flux', 
          'kplus_Flux', 
          'kzero_Flux', 
          'nucleoninexsec_Flux', 
          'nucleonqexsec_Flux', 
          'nucleontotxsec_Flux', 
          'piminus_Flux', 
          'pioninexsec_Flux', 
          'pionqexsec_Flux', 
          'piontotxsec_Flux', 
          'piplus_Flux', 
          'total_multisim', 
          'total_reweighted', 
          'total_flux']

labels = ['Horn Optics',
          r'$\pi^+$ Flux',
          r'$K^+$ Flux',
          'Other',
          'Total multisim',
          'Total',
          'Total Flux',]

for i,state_name in enumerate(state_names):
  hist = uproot.open(f'{state_folder}/{state_fnames[i]}')
  #GEt bin values and edges
  keys = hist.keys()
  edges = hist[keys[0]].axis().edges() #They all have the same bin widths
  values = np.zeros((len(keys),len(edges)-1)) #One less bin edge shape
  for j,key in enumerate(keys): #Iterate over all universes
    #print(hist[key].values()[5],',')
    values[j] = hist[key].values()
  #values*=pot_normalization
  if state_fnames[i] == 'state_Flux_multisim.root':
    nom = values[0]
  #Calculate systematic uncertainties with bin heights. We can use this to make an error band

  (height_constrained,height,fracunc_constrained,fracunc,
  fracbias_constrained,fracbias,fracstat_constrained,fracstat) = helpers.frac_syst_stat_bias_heights(values,edges,W_arr_stat)
  
  #Store in arrays as well, we probably won't use the dict
  heights_constrained.append(height_constrained)
  heights.append(height)
  fracuncs_constrained.append(fracunc_constrained)
  fracuncs.append(fracunc)
  fracbiases_constrained.append(fracbias_constrained)
  fracbiases.append(fracbias)
  fracstats_constrained.append(fracstat_constrained)
  fracstats.append(fracstat)

  #print(f'{state_name} : ',np.diff(values).shape,np.diff(values),edges)
#***----Custom setting!---***
#edges[-1] = 3
print(edges)
#print(np.shape(heights_constrained),np.shape(height_constrained))
#Optics
(opticsheights_constrained,opticsheights,opticsuncs_constrained,opticsuncs,opticsbiases_constrained,
 opticsbiases,opticsstats_constrained,opticsstats) = combine_frac_height_vals(optics_indeces,heights,heights_constrained,
                      fracuncs_constrained,fracuncs,fracbiases_constrained,
                      fracbiases,fracstats_constrained,fracstats,return_combined=True)

#Piplus
(piplusheights_constrained,piplusheights,piplusuncs_constrained,piplusuncs,piplusbiases_constrained,
 piplusbiases,piplusstats_constrained,piplusstats) = combine_frac_height_vals(piplus_index,heights,heights_constrained,
                      fracuncs_constrained,fracuncs,fracbiases_constrained,
                      fracbiases,fracstats_constrained,fracstats,return_combined=True)
#Kplus
(kplusheights_constrained,kplusheights,kplusuncs_constrained,kplusuncs,kplusbiases_constrained,
 kplusbiases,kplusstats_constrained,kplusstats) = combine_frac_height_vals(kplus_index,heights,heights_constrained,
                      fracuncs_constrained,fracuncs,fracbiases_constrained,
                      fracbiases,fracstats_constrained,fracstats,return_combined=True)
#Other
(otherheights_constrained,otherheights,otheruncs_constrained,otheruncs,otherbiases_constrained,
 otherbiases,otherstats_constrained,otherstats) = combine_frac_height_vals(other_indeces,heights,heights_constrained,
                      fracuncs_constrained,fracuncs,fracbiases_constrained,
                      fracbiases,fracstats_constrained,fracstats,return_combined=True)

#Total
(fluxheights_constrained,fluxheights,fluxuncs_constrained,fluxuncs,fluxbiases_constrained,
 fluxbiases,fluxstats_constrained,fluxstats) = combine_frac_height_vals(flux_index,heights,heights_constrained,
                      fracuncs_constrained,fracuncs,fracbiases_constrained,
                      fracbiases,fracstats_constrained,fracstats,return_combined=True)

#Combine into lists
labels = ['Horn optics',r'$\pi^+$ flux',r'$K^+$ flux','Other','Total']
heights = np.array([opticsheights,piplusheights,kplusheights,otherheights,fluxheights])
fracuncs = np.array([opticsuncs,piplusuncs,kplusuncs,otheruncs,fluxuncs])
fracbiases = np.array([opticsbiases,piplusbiases,kplusbiases,otherbiases,fluxbiases])
colors = cm.get_cmap('viridis', len(labels)-1)

fig,ax,ax_unc = plot_resolution_bias_stat(nom,heights,edges,labels,colors,fracunc=fracuncs,fracbias=None,fracunc_constrained=None,
                                        fracbias_constrained=None,fracstat=None,fracstat_constrained=None)
#ax.set_ylim([0,None])
#ax_unc.set_ylim([0,0.2])
ax.set_xlim([0,3])
ax_unc.set_xlim([0,3])

ax.set_title('Flux Uncertainties')
ax_unc.set_xlabel(r'True $E_\nu$ (GeV)')
plotters.save_plot('flux_syst',folder_name=plot_folder,fig=fig)
ax.set_yscale('log')
plotters.save_plot('flux_syst_logy',folder_name=plot_folder,fig=fig)

fig,ax,ax_unc = plot_resolution_bias_stat(nom,heights,edges,labels,colors,fracunc=fracuncs,fracbias=fracbiases,fracunc_constrained=None,
                                        fracbias_constrained=None,fracstat=None,fracstat_constrained=None)
ax.set_xlim([0,3])
ax_unc.set_xlim([0,3])
ax.set_title('Flux Uncertainties w/ Bias')
ax_unc.set_xlabel(r'True $E_\nu$ (GeV)')
plotters.save_plot('flux_syst_wbias',folder_name=plot_folder,fig=fig)
ax.set_yscale('log')
plotters.save_plot('flux_syst_wbias_logy',folder_name=plot_folder,fig=fig)

#print([fname[6:-5] for fname in fnames if 'band' not in fname and 'nom' not in fname])

#weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_systs_truth_cuts'
#weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_3_systs_nuescat_cut'
#weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_29_systs_truth_cuts_eleeng'
weights_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_31_systs_truth_cuts_eleeng'
#'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_30_systs_5bins'

#sample_pot = 1.55131e+19 #Size of sample used to make
#sample_pot = 6.16911e+19
sample_pot = 3.08e20 #Total from MCP2022A
nom_pot = 10e21 #Expected POT of sbnd
pot_sample_size = 1e9 #Renormalize x-axis to this value for nice numbers

pot_normalization = 1/sample_pot #Get normalization to reach nominal flux



