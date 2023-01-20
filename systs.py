import numpy as np
import uproot
import pandas as pd
import matplotlib.pyplot as plt
import helpers
plt.style.use(['science','no-latex'])

sample = 'CCNuE'

day = helpers.day

genie_syst_tree = uproot.open(f'/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/2023_1_16_all_lowpot/2023_1_16/error_band_interaction_{sample}.root:systs;1')
genie_syst = genie_syst_tree.arrays(library='pd')

flux_syst_tree = uproot.open(f'/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/2023_1_16_all_lowpot/2023_1_16/error_band_flux_{sample}.root:systs;1')
flux_syst = flux_syst_tree.arrays(library='pd')

genie_syst.loc[:,'yerr_minus'] = (genie_syst.loc[:,'ynom']-genie_syst.loc[:,'y0'])/genie_syst.loc[:,'ynom']
genie_syst.loc[:,'yerr_plus'] = abs((genie_syst.loc[:,'ynom']-genie_syst.loc[:,'y1'])/genie_syst.loc[:,'ynom'])

flux_syst.loc[:,'yerr_minus'] = (flux_syst.loc[:,'ynom']-flux_syst.loc[:,'y0'])/flux_syst.loc[:,'ynom']
flux_syst.loc[:,'yerr_plus'] = abs((flux_syst.loc[:,'ynom']-flux_syst.loc[:,'y1'])/flux_syst.loc[:,'ynom'])

genie_syst = genie_syst.fillna(0)
flux_syst = flux_syst.fillna(0)

x = genie_syst.loc[:,'xnom'].values

y_label = 'yerr_minus'
y = genie_syst.loc[:,y_label].values
y2 = flux_syst.loc[:,y_label].values

y1_label = 'yerr_plus'
y1 = genie_syst.loc[:,y1_label].values
y3 = flux_syst.loc[:,y1_label].values
fig, ax = plt.subplots(figsize=(8, 2), tight_layout=True)

flux_avg = np.mean((y2+y3)/2) #Change if bin sizes are different
genie_avg = np.mean((y+y1)/2) #Change if bin sizes are different

ax.plot(x,(y+y1)/2,label=f'genie ({genie_avg*100:.1f}%)',color='red')
ax.plot(x,(y2+y3)/2,label=f'flux ({flux_avg*100:.1f}%)',color='blue')
ax.legend(fontsize=16)
ax.set_xlabel('Energy',fontsize=14)
ax.set_ylabel('Fractional Uncertainty',fontsize=14)
#ax.set_title('Fractional Systematic Uncertainties',fontsize=18)

#plt.savefig('systs.png')

#Error bar hist
x = genie_syst.loc[:,'xnom']
dx = genie_syst.loc[:,'dx'].values[0]

y_syst = genie_syst.loc[:,'ynom']
y_flux = flux_syst.loc[:,'ynom']

y0_syst = genie_syst.loc[:,'y0']
y0_flux = flux_syst.loc[:,'y0']

y1_syst = genie_syst.loc[:,'y1']
y1_flux = flux_syst.loc[:,'y1']

yerrs_syst = [y1_syst,y0_syst]
yerrs_flux = [y1_flux,y0_flux]

helpers.save_plot(f'frac_systs_{sample}',folder_name=f'/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/2023_1_16_all_lowpot/systs')

fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)

ax.bar(x,y_flux,width=dx,color='blue',alpha=0.3,edgecolor='blue',linewidth=2)
ax.bar(x,height=y1_flux-y0_flux,width=dx,bottom=y0_flux,alpha=0.4)
ax.set_ylim([0,None])
ax.set_xlabel('Energy [GeV]',fontsize=16)
ax.set_ylabel('Counts',fontsize=16)
helpers.set_style(ax)
ax.set_title('Flux Uncertainties',fontsize=20)
helpers.save_plot(f'flux_systs_{sample}',folder_name=f'/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/2023_1_16_all_lowpot/systs')

fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)

ax.bar(x,y_syst,width=dx,color='red',alpha=0.3,edgecolor='red',linewidth=2)
ax.bar(x,height=y1_syst-y0_syst,width=dx,bottom=y0_syst,alpha=0.4,color='red')
ax.set_ylim([0,None])
ax.set_xlabel('Energy [GeV]',fontsize=16)
ax.set_ylabel('Counts',fontsize=16)
helpers.set_style(ax)
ax.set_title('Interaction Uncertainties',fontsize=20)
helpers.save_plot(f'genie_systs_{sample}',folder_name=f'/sbnd/data/users/brindenc/analyze_sbnd/nue/plots/2022A/2023_1_16_all_lowpot/systs')


