import pandas as pd
import seaborn as sns
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time
import matplotlib
matplotlib.use('Agg')

import sys
sys.path.append('../')
from CAFdata import *
import helpers
from datetime import date

day = date.today().strftime("%Y_%m_%d")
NuEScat_dir = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_27_ethetacut'
surname = ''
fnames = [f'cut_events_basic{surname}.root',f'cut_events_trk{surname}.root',f'cut_events_shw{surname}.root']
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

def plot_all_keys(df,ftype='.png'):
  """
  plot_all_keys(df)

  Parameters:
  df (pandas dataframe): dataframe that contains the data to be plotted

  Returns:
  None
"""
  counter = 0
  folder_name=f'{NuEScat_dir}/plots'
  fnames = []
  alphas = [1,0.4,0.4,0.4,0.4,0.4]
  os.system(f'mkdir -p {folder_name}')
  for key_x in df.keys():
    for key_y in df.keys():
      if key_x != key_y or key_x == 'evt_type': continue
      xs = [] #x axis points
      ys = [] #y axis points
      labels = []
      alphas = np.zeros(len(event_names))
      fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)
      for i,name in enumerate(event_names):
        if name == r'$\nu + e$':
          alphas[i] = 1
        else: 
          alphas[i] = 0.01
        xs.append(df.loc[df['evt_type']==name,key_x].values)
        ys.append(df.loc[df['evt_type']==name,key_y].values)
        labels.append(f'{name} ({len(xs[i])})')
        if key_x != key_y:
          sns.scatterplot(x=xs[i],y=ys[i],label=labels[i],ax=ax,alpha=alphas[i])
          ax.set_xlabel(f'{key_x}')
          ax.set_ylabel(f'{key_y}')
        if key_x == key_y and key_x != 'evt_type':
          #weights = [abs(x/len(x)) for x in xs] #Weight each bin to its weight. Each histogram bin should be one
          #ax.hist(weights,label=labels,alpha=0.3,edgecolor=colors[i],stacked=True,density=False)#,weights=weights)
          #sns.histplot(x=xs,ax=ax,label=labels,alpha=alphas,common_norm=False,
          #common_bins=False,bins=10,color=colors)
          lm = sns.kdeplot(x=xs[i],ax=ax,label=labels[i],warn_singular=False,alpha=alphas[i])
          #ax.set_xlabel(f'{key_x}')
          #ax.set_ylabel('Counts')
        
        
      ax.legend()
      helpers.set_style(ax,legend_size=12)
      fnames.append(f'{key_x.replace(".","_")}_{key_y.replace(".","_")}_dist')
      plt.savefig(f'{fnames[counter]}{ftype}',bbox_inches = "tight",dpi=300)
      print(fnames[counter])
      plt.clf()
      plt.close('all')
      del ax
      del fig
      counter+=1
  for fname in fnames:
    os.system("mv " + fname + ftype + f' {folder_name}/')

def pandora_reco_cut(events):
  """
  keep following pairs of reco objects (nshw,ntrk,stub):
  - 1,1,0
  - 1,0,1
  - 2,0,0
  - 1,0,0
  - 0,0,1
  - 0,1,0
  - 0,1,1
  """
  return events

#Use uproot to load all the ttrees and concatenate them
dfs = []
for fname in fnames:
  #print(fname)
  tree = uproot.open(f'{NuEScat_dir}/{fname}:rectree;1')
  df = helpers.get_df(tree,tree.keys(),hdrkeys=['run','subrun','evt'])
  dfs.append(df)
events = pd.concat(dfs,axis=1)

#Mask event type and get their names
events = set_event_type(events)
event_names = categories

#Isolate data from different reconstructed objects, shws trks
lshw_keys = [key for key in events.keys() if key[:5] == 'lshw.']
lshw_keys.extend(['evt_type'])
slshw_keys = [key for key in events.keys() if key[:6] == 'slshw.']
slshw_keys.extend(['evt_type'])

lshw = events.loc[:,lshw_keys]
slshw = events.loc[:,slshw_keys]

lshw = lshw[(lshw != -9999).all(1)] #Remove all -999 dummy values
lshw = lshw[(lshw != -999).all(1)] #Remove all -999 dummy values
lshw = lshw[(lshw != -5).all(1)] #Remove all -5 dummy values
lshw = lshw.loc[lshw.loc[:,'lshw.dedx']<100] #Remove all huge dedx values

slshw = slshw[(slshw != -9999).all(1)] #Remove all -999 dummy values
slshw = slshw[(slshw != -999).all(1)] #Remove all -999 dummy values
slshw = slshw[(slshw != -5).all(1)] #Remove all -5 dummy values
slshw = slshw.loc[slshw.loc[:,'slshw.dedx']<100] #Remove all huge dedx values

ltrk_keys = [key for key in events.keys() if key[:5] == 'ltrk.']
ltrk_keys.extend(['evt_type'])
sltrk_keys = [key for key in events.keys() if key[:6] == 'sltrk.']
sltrk_keys.extend(['evt_type'])

ltrk = events.loc[:,ltrk_keys]
sltrk = events.loc[:,sltrk_keys]

ltrk = ltrk[(ltrk != -9999).all(1)] #Remove all -999 dummy values
ltrk = ltrk[(ltrk != -999).all(1)] #Remove all -999 dummy values
ltrk = ltrk[(ltrk != -5).all(1)] #Remove all -5 dummy values
#ltrk = ltrk.loc[ltrk.loc[:,'ltrk.dedx']<100] #Remove all huge dedx values

sltrk = sltrk[(sltrk != -9999).all(1)] #Remove all -999 dummy values
sltrk = sltrk[(sltrk != -999).all(1)] #Remove all -999 dummy values
sltrk = sltrk[(sltrk != -5).all(1)] #Remove all -5 dummy values
#sltrk = sltrk.loc[sltrk.loc[:,'sltrk.dedx']<100] #Remove all huge dedx values

other_keys = []
for key in events.keys():
  if key not in slshw_keys and key not in sltrk_keys and key not in ltrk_keys and key not in lshw_keys:
    if key != 'run' and key != 'subrun' and key != 'evt':
      other_keys.append(key)
other_keys.extend(['evt_type'])
other_events = events.loc[:,other_keys]
other_events = other_events[(other_events != -9999).all(1)] #Remove all -999 dummy values
other_events = other_events[(other_events != 9999).all(1)] #Remove all -999 dummy values
other_events = other_events[(other_events != -999).all(1)] #Remove all -999 dummy values


#print(other_keys)
plot_all_keys(other_events)
plot_all_keys(lshw)
plot_all_keys(slshw)
plot_all_keys(ltrk)
plot_all_keys(sltrk)

