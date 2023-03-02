import pandas as pd
import seaborn as sns
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

import sys
sys.path.append('../')
sys.path.append('/sbnd/app/users/brindenc/mypython/bc_utils')
from utils import plotters,pic
from CAFdata import *
import helpers
import cuts
from datetime import date

day = date.today().strftime("%Y_%m_%d")
plotters.use_science_style()

NuEScat_dir = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_1_fullsample_few_recocut_recoslc'
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


# def select_electron(events,objs=['lshw','ltrk','slshw','sltrk']):
#   """
#   Function to select which reco object (lshw,ltrk,slshw,etc.) is the electron
#   We will use the angle as our main cut, but also ensure that the number of points is fine
  
#   Parameters:
#   events (dataframe): Contains event information 
#   """
#   for ind,row in events.iterrows():



#Use uproot to load all the ttrees and concatenate them
dfs = []
for fname in fnames:
  #print(fname)
  tree = uproot.open(f'{NuEScat_dir}/{fname}:rectree;1')
  print(fname)
  df = helpers.get_df(tree,tree.keys(),hdrkeys=['run','subrun','evt'],library='np')
  print(df['run'],df['Etheta'])
  df = pd.DataFrame.from_dict(df)
  #df.set_index(['run','subrun','evt'])
  #df.sort_index()
  #dfs.append(df)
events = pd.concat(dfs,axis=1)
events.iloc[0].to_csv('head.csv')

#Mask event type and get their names
events = set_event_type(events)
events.loc[:,'nreco'] = events.nshw + events.ntrk

events = cuts.apply_cuts(events)

print('Trying to save')
events.to_pickle(f'{NuEScat_dir}/events_raw.pkl')


