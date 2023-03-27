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

#NuEScat_dir = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_1_fullsample_few_recocut_recoslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_8_fullsample_nue_few_recocut_recoslc_theta'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_14_pure_nue_truefv_recoslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_14_fullsample_softetheta_recoslc'
NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_15_pure_nue_truefv_recoslc'
surname = ''
fnames = [f'cut_events_basic{surname}.root',f'cut_events_trk{surname}.root',f'cut_events_shw{surname}.root']



#Use uproot to load all the ttrees and concatenate them
dfs = []
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
#events.iloc[0].to_csv('head.csv')

#Mask event type and get their names
events = helpers.set_event_type(events)
events.loc[:,'nreco'] = events.nshw + events.ntrk

#Get electron


#events = cuts.apply_cuts(events)

print('Trying to save')
events.to_csv(f'{NuEScat_dir}/events_raw.csv')


