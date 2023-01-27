import pandas as pd
import seaborn as sns
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import sys
sys.path.append('../')
from CAFdata import *
import helpers
from datetime import date

day = date.today().strftime("%Y_%m_%d")
NuEScat_dir = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_26'
fnames = ['cut_events_basic_stride10.root','cut_events_trk_stride10.root','cut_events_shw_stride10.root']

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
  events.loc[events.loc[:,'evt_type'] == 0,'evt_type'] = 'NuEScat'
  events.loc[events.loc[:,'evt_type'] == 1,'evt_type'] = 'NCPi0'
  events.loc[events.loc[:,'evt_type'] == 2,'evt_type'] = 'NC'
  events.loc[events.loc[:,'evt_type'] == 3,'evt_type'] ='CCNuMu'
  events.loc[events.loc[:,'evt_type'] == 4,'evt_type'] = 'CCNuE'
  events.loc[events.loc[:,'evt_type'] == 5,'evt_type'] = 'Dirt'
  events.loc[events.loc[:,'evt_type'] == 6,'evt_type'] ='Other'
  return events

def plot_all_keys(df):
  """
  plot_all_keys(df)

  This function takes a dataframe (df) as input and generates a scatter plot 
  or a kde plot for all possible combinations of keys from the dataframe. 
  The function loops through all keys of the dataframe and for each pair of 
  keys, it creates a figure and axes, set the style of the plot, separates 
  the data by event name, and plots the data. The plot can be either a 
  scatter plot or a kde plot, depending on whether the keys are the same 
  or not. Each plot is saved to a specified folder with the keys as part 
  of the file name.

  Parameters:
  df (pandas dataframe): dataframe that contains the data to be plotted

  Returns:
  None
"""
  for key_x in df.keys():
    for key_y in df.keys():
      fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)
      helpers.set_style(ax)
      for name in event_names:
        x = df.loc[df['evt_type']==name,key_x]
        y = df.loc[df['evt_type']==name,key_y]
        if key_x == key_y and key_x!='evt_type' and key_y!= 'evt_type':
          sns.kdeplot(x=x,ax=ax,label=f'{name} ({len(x)})')
          #ax.hist(x,label=f'{name} ({len(x)})',fill=None)
        else:
          sns.scatterplot(x=x,y=y,label=f'{name} ({len(x)})',ax=ax,alpha=0.5)
      print(f'{key_x} {key_y}')
      ax.legend();
      helpers.save_plot(f'{key_x}_{key_y}',folder_name=f'{NuEScat_dir}/plots')
      plt.close(fig)
      del ax
      del fig

#Use uproot to load all the ttrees and concatenate them
dfs = []
for fname in fnames:
  tree = uproot.open(f'{NuEScat_dir}/{fname}:rectree;1')
  df = helpers.get_df(tree,tree.keys(),hdrkeys=['run','subrun','evt'])
  dfs.append(df)
events = pd.concat(dfs,axis=1)

#Mask event type and get their names
events = set_event_type(events)
event_names = events['evt_type'].unique()
colors = cm.get_cmap('Spectral', len(event_names))

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


print(other_keys)
plot_all_keys(other_events)
#plot_all_keys(lshw)
#plot_all_keys(slshw)
#plot_all_keys(ltrk)
#plot_all_keys(sltrk)

