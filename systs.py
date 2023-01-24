import numpy as np
import uproot
import pandas as pd
import matplotlib.pyplot as plt
import helpers
import os
import matplotlib.cm as cm
plt.style.use(['science','no-latex'])

sample = ''

day = helpers.day

state_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_23_test'
plot_folder = '/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/2023_1_23_test'

fnames = os.listdir(state_folder)
colors = cm.get_cmap('Spectral', len(fnames))
trees = []
syst_dfs = []
labels = []
for i,fname in enumerate(fnames):
  tree = uproot.open(f'{state_folder}/{fname}:systs;1')
  df = tree.arrays(library='pd')
  df.loc[:,'yerr_minus'] = (df.loc[:,'ynom']-df.loc[:,'y0'])/df.loc[:,'ynom']
  df.loc[:,'yerr_plus'] = (df.loc[:,'y1']-df.loc[:,'ynom'])/df.loc[:,'ynom']
  df = df.fillna(0)
  labels.append(fname[11:-5])
  syst_dfs.append(df)

fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)
fig2, ax2 = plt.subplots(figsize=(8, 2), tight_layout=True)
helpers.set_style(ax)
helpers.set_style(ax2)

for i,df in enumerate(syst_dfs):
  fig3, ax3 = plt.subplots(figsize=(8, 2), tight_layout=True)
  helpers.set_style(ax3)
  x = df.loc[:,'xnom'].values
  dx = df.loc[:,'dx'].values

  y0_err = df.loc[:,'yerr_minus'].values
  y1_err = df.loc[:,'yerr_plus'].values
  y = df.loc[:,'ynom'].values
  y1 = df.loc[:,'y1'].values
  y0 = df.loc[:,'y0'].values

  avg = np.mean((y0_err+y1_err)/2)
  if labels[i] == 'GENIE_multisim':
    lw = 2
    ax2.bar(x,(y0_err+y1_err)/2,width=dx,label=f'{labels[i]} ({avg*100:.1f}%)',linewidth=lw,
    fill=False,edgecolor=colors(i))
  elif labels[i] == 'Flux_multisim':
    lw = 2
    ax2.bar(x,(y0_err+y1_err)/2,width=dx,label=f'{labels[i]} ({avg*100:.1f}%)',linewidth=lw,
    fill=False,edgecolor=colors(i))
  ax.bar(x,(y0_err+y1_err)/2,width=dx,label=f'{labels[i]} ({avg*100:.1f}%)',linewidth=lw,
    fill=False,edgecolor=colors(i))
  ax3.bar(x,y,width=dx,color=colors(i),fill=False,edgecolor=colors(i),linewidth=2,
  label=f'{labels[i]}')
  ax3.bar(x,height=y1-y0,width=dx,bottom=y0,alpha=0.4,color=colors(i))
  #ax3.set_title('Interaction Uncertainties',fontsize=20)
  ax3.legend(fontsize=14)
  ax3.set_xlabel('Energy [GeV]',fontsize=14)
  ax3.set_ylabel('Counts',fontsize=14)
  helpers.save_plot(f'band_{labels[i]}',folder_name=f'{plot_folder}/systs')

# ax.plot(x,(y+y1)/2,label=f'genie ({genie_avg*100:.1f}%)',color='red')
# ax.plot(x,(y2+y3)/2,label=f'flux ({flux_avg*100:.1f}%)',color='blue')
ax.legend(fontsize=12)
ax2.legend(fontsize=12)
ax.set_xlabel('Energy',fontsize=16)
ax2.set_xlabel('Energy',fontsize=16)
ax.set_ylabel('Fractional Uncertainty',fontsize=16)
ax2.set_ylabel('Fractional Uncertainty',fontsize=16)
helpers.save_plot(f'frac_systs_all',folder_name=f'{plot_folder}/systs',fig=fig)
helpers.save_plot(f'frac_systs',folder_name=f'{plot_folder}/systs',fig=fig2)
#ax.set_title('Fractional Systematic Uncertainties',fontsize=18)



