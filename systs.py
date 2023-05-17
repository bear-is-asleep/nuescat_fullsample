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

#suffix = '_systs_nuescat_cut'
#state_folder = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_3_systs_nuescat_cut'
#state_folder = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_systs_truth_cuts'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_16_reweight_flux'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_29_systs_truth_cuts_eleeng'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_29_reweight_flux_stride100'
#state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_30_reweight_flux_var'
state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_31_systs_truth_cuts_eleeng'

plot_folder = f'{state_folder}/plots'

fnames = os.listdir(state_folder)
colors = cm.get_cmap('Spectral', len(fnames))
trees = []
syst_dfs = []
labels = []
for i,fname in enumerate(fnames):
  if fname[:5] == 'error':
    tree = uproot.open(f'{state_folder}/{fname}:systs;1')
    df = tree.arrays(library='pd')
    df.loc[:,'yerr_minus'] = (df.loc[:,'ynom']-df.loc[:,'y0'])/df.loc[:,'ynom']
    df.loc[:,'yerr_plus'] = (df.loc[:,'y1']-df.loc[:,'ynom'])/df.loc[:,'ynom']
    df = df.fillna(0)
    labels.append(fname[11:-5])
    syst_dfs.append(df)

fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)
fig2, ax2 = plt.subplots(figsize=(8, 4), tight_layout=True)
fig4,ax4 = plt.subplots(figsize=(8,6),tight_layout=True)
ax.set_xlim([None,8])
ax2.set_xlim([None,8])

rects = [] #Stores bar chart labels for horizontal error plot

for i,df in enumerate(syst_dfs):
  fig3, ax3 = plt.subplots(figsize=(8, 4), tight_layout=True)
  ax3.set_xlim([None,8])
  x = df.loc[:,'xnom'].values[1:-1]
  dx = df.loc[:,'dx'].values[1:-1]

  y0_err = df.loc[:,'yerr_minus'].values[1:-1]
  y1_err = df.loc[:,'yerr_plus'].values[1:-1]
  

  y = df.loc[:,'ynom'].values[1:-1]
  y1 = df.loc[:,'y1'].values[1:-1]
  y0 = df.loc[:,'y0'].values[1:-1]

  avg = np.mean((y0_err+y1_err)/2)
  lw = 2
  #These uncertainties are the main ones, interaction, flux, total
  if labels[i] == 'GENIE_multisim':
    lw = 2
    ax2.bar(x,(y0_err+y1_err)/2,width=dx,label=f'{labels[i]} ({avg*100:.1f}%)',linewidth=lw,
    fill=False,edgecolor=colors(0))
  elif labels[i] == 'Flux_multisim':
    lw = 2
    ax2.bar(x,(y0_err+y1_err)/2,width=dx,label=f'{labels[i]} ({avg*100:.1f}%)',linewidth=lw,
    fill=False,edgecolor=colors(5))
  elif labels[i] == 'total_multisim':
    #print(y1_err,y0_err,(y0_err+y1_err)/2)
    lw = 2
    #ax2.step(x,(y0_err+y1_err)/2,label=f'{labels[i]} ({avg*100:.1f}%)',linewidth=lw,where='post',color=colors(10))
    ax2.bar(x,(y0_err+y1_err)/2,width=dx,label=f'{labels[i]} ({avg*100:.1f}%)',linewidth=lw,
    fill=False,edgecolor=colors(10))
  #Add bar for all syst uncertainties
  ax.bar(x,(y0_err+y1_err)/2,width=dx,label=f'{labels[i]} ({avg*100:.1f}%)',linewidth=lw,
    fill=False,edgecolor=colors(i))
  #Individual systematic uncertainty
  ax3.bar(x,
    y,
    width=dx,
    color=colors(i),
    fill=False,
    edgecolor=colors(i),
    linewidth=2,
    label=f'{labels[i]}'
  )
  ax3.bar(x,
    height=y1-y0,
    width=dx,
    bottom=y0,
    alpha=0.4,
    color=colors(i)
  )
  plotters.set_style(ax3,legend_size=16,axis_size=14)
  #ax3.set_title('Interaction Uncertainties',fontsize=20)
  ax3.legend()
  ax3.set_xlabel('Energy [GeV]')
  ax3.set_ylabel('Counts')
  plotters.save_plot(f'band_{labels[i]}',folder_name=f'{plot_folder}/systs')

  #Horizontal bar chart
  rect = ax4.bar(f'{labels[i]}',# ({avg*100:.1f}%)',
    height=avg*200/2, #times two to encompass +- values
    width=0.6,
    bottom=-np.mean(y0_err)*100/2,
    color=colors(i)
  )

# ax.plot(x,(y+y1)/2,label=f'genie ({genie_avg*100:.1f}%)',color='red')
# ax.plot(x,(y2+y3)/2,label=f'flux ({flux_avg*100:.1f}%)',color='blue')
ax.legend()
ax.set_xlabel('True energy [GeV]')
ax.set_ylabel('Fractional Uncertainty')
plotters.set_style(ax,legend_size=12,axis_size=16)
plotters.save_plot(f'frac_systs_all',folder_name=f'{plot_folder}/systs',fig=fig)

ax2.legend()
ax2.set_xlabel('Energy')
ax2.set_ylabel('Fractional Uncertainty')
plotters.set_style(ax2,legend_size=12,axis_size=16)
plotters.save_plot(f'frac_systs',folder_name=f'{plot_folder}/systs',fig=fig2)

#ax4.legend()
ax4.set_ylabel('Fractional Uncertainty %')
fig4.autofmt_xdate(rotation=60)
ax4.grid(True)
plotters.set_style(ax4,axis_size=16)
plotters.save_plot(f'horizontal_systs',folder_name=f'{plot_folder}/systs',fig=fig4)
#ax.set_title('Fractional Systematic Uncertainties',fontsize=18)



