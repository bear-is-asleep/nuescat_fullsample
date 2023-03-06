import sys
sys.path.append('../')
sys.path.append('/sbnd/app/users/brindenc/mypython/bc_utils')
import cuts
import pandas as pd
import seaborn as sns
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from time import time

from utils import plotters,pic
from CAFdata import *
import helpers
from datetime import date

day = date.today().strftime("%Y_%m_%d")
plotters.use_science_style()
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_2_27_fullsample_nue_no_cuts_truthslc'
NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_1_fullsample_few_recocut_recoslc'
surname = ''
fnames = [f'cut_events_basic{surname}.root',f'cut_events_trk{surname}.root',f'cut_events_shw{surname}.root']
folder_name=f'{NuEScat_dir}/plots'

colors = ['blue','purple','yellow','red','green','brown','black']
categories = [r'$\nu + e$',r'NC $\pi^0$','NC',r'CC$\nu_\mu$',r'CC$\nu_e$','Dirt','Other']
def set_event_type(events):
  """
  set_event_type(events)

  This function takes a dataframe (events) as input, and sets the 
  'evt_type' column to a string value based on the numerical value of
  'evt_type' column. The numerical values are mapped to string values: 
  0: 'NuEScat', 1: 'NCPi0', 2: 'NC', 3: 'CCNuMu', 4: 'CCNuE', 5: 'Dirt', 6: 'Other'.

  Parameters:
  events (pandas dataframe): dataframe that contains the event data

  Returns:
  events (pandas dataframe): modified dataframe with the 'evt_type' column as string values
  """
  
  for i,cat in enumerate(categories):
    events.loc[events.loc[:,'evt_type'] == i,'evt_type'] = cat
  return events

#Use uproot to load all the ttrees and concatenate them
dfs = []
for fname in fnames:
  print(fname)
  tree = uproot.open(f'{NuEScat_dir}/{fname}:rectree;1')
  #use numpy library to convert numpy array to dataframe
  df = helpers.get_df(tree,
                      tree.keys(),
                      hdrkeys=['run','subrun','evt'],
                      library='np')
  dfs.append(df)
events = pd.concat(dfs,axis=1)

print('-----Cleaning data')
#Mask event type and get their names
events = set_event_type(events)
events.loc[:,'nreco'] = events.nshw + events.ntrk
#events = cuts.apply_cuts(events)
event_names = categories

print('-----Make efficiency plot')
#Cut values
cut_ranges = [np.arange(0.0001,0.02,0.0001),
              np.arange(10,400,10)]

#Cut functions
cut_funcs = [cuts.cEtheta2,
             cuts.cLShwLen]

#Cut keys
cut_keys = ['Etheta',
            'lshw.len']

#Bool keys - says whether event passes or fails cut
bool_keys = ['cEtheta2',
             'cLShwLen']

for j,cut_func in enumerate(cut_funcs):
  print(cut_keys[j])
  cut_vals = cut_ranges[j]
  #Efficiency,purity and f1 initialization
  effs = np.zeros(cut_vals.shape[0])
  purs = effs.copy()
  f1s = effs.copy()
  for i,cut_val in enumerate(cut_vals):
    #Apply cut from list of cut functions
    events = cut_func(events,
                      cut_key=cut_keys[j],
                      bool_key=bool_keys[j],
                      cut_val=cut_val)
    #print(events.head())
    mask = events.loc[:,bool_keys[j]].values #get events that pass cut
    pur,eff,f1 = cuts.efficiency_purity_f1_calc(events,mask)
    effs[i] = eff
    purs[i] = pur
    f1s[i] = f1

  fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)
  ax2 = ax.twinx()
  ln1 = ax.plot(cut_vals,effs,label='Eff')
  ln2 = ax.plot(cut_vals,purs,label='Pur',color='red')
  ln3 = ax2.plot(cut_vals,f1s,label='F1',color='green')

  ax.set_xlabel(f'{cut_keys[j]} cut value')
  ax.set_ylabel('Pur and Eff')
  ax2.set_ylabel(r'F1 = $\frac{1}{2} (1/$pur$ + 1/ $eff$ )$')
  ax.set_title(f'')

  plotters.set_style(ax)
  plotters.set_style(ax2)

  #combine legends
  lines1, labels1 = ax.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  lines = lines1 + lines2
  labels = labels1 + labels2
  ax.legend(lines, labels)
  plotters.save_plot(f'eff_pur_f1_{bool_keys[j]}',folder_name=folder_name+'/effs')






