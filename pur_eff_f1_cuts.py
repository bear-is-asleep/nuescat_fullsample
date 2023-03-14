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

import logging
logging.basicConfig(level=logging.INFO)

make_plots = True
apply_cuts = True

day = date.today().strftime("%Y_%m_%d")
plotters.use_science_style()
NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_2_27_fullsample_nue_no_cuts_truthslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_1_fullsample_few_recocut_recoslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_9_fullsample_few_recocut_recoslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_13_fullsample_nue_few_recocut_recoslc'
os.system(f'mkdir -p {NuEScat_dir}/plots')
os.system(f'mkdir -p {NuEScat_dir}/plots/effs')

surname = ''
fnames = [f'cut_events_basic{surname}.root',f'cut_events_trk{surname}.root',f'cut_events_shw{surname}.root']
folder_name=f'{NuEScat_dir}/plots'

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
events = helpers.set_event_type(events)
events.loc[:,'nreco'] = events.nshw + events.ntrk
#events = cuts.mask_cosmics(events,time_threshold=[0.,1.6]) #Mask all events that have crt activity
events = cuts.select_electron(events) #Should just be lshw

print('-----Make efficiency plots')
#Cut values
etheta_range = np.arange(0.0001,0.01,step=0.001)
trk_len_range = np.arange(2,400,step=2)
trk_eng_range = np.arange(0.01,1,step=0.01)
openangle_range = np.arange(1e-5,45,step=1)
bdt_score_range = np.arange(0,1.05,step=0.05)
dedx_range = np.arange(1,10,step=0.1)

"""
Selection dictionary has the cut key as the key and in it is
cut_range : the range of cut values to test
cut_function : what type of cut
cut_key : which key it is in the df
equality : cut on values lessthan, greaterthan, equal?

Next 3 are for NReco cut only, set to None if using a box cut
nreco_key : which key to check for nreco objects, i.e. nshw
nreco : number of reco objects, i.e. nshw = 1 would be 1
nreco_equality : check for reco objects lessthan, greaterthan, equal, i.e. nshw = 1 would be equal
"""



nue_sel_dict = {'cLTrkLen': {
                'cut_range':trk_len_range,
                'cut_function':cuts.cBoxCut,
                'cut_key':'ltrk.len',
                'equality':'lessthan',
                'nreco_key':None,
                'nreco':None,
                'nreco_equality':None
                },
                'cLTrkEng': {
                'cut_range':trk_eng_range,
                'cut_function':cuts.cBoxCut,
                'cut_key':'ltrk.eng',
                'equality':'lessthan',
                'nreco_key':None,
                'nreco':None,
                'nreco_equality':None
                },
}

# cut_ranges = [etheta_range,
#               len_range,
#               openangle_range,
#               bdt_score_range,
#               dedx_range,
#               etheta_range,
#               len_range,
#               len_range]

# #Cut functions
# cut_funcs = [cuts.cBoxCut,
#              cuts.cBoxCut,
#              cuts.cBoxCut,
#              cuts.cBoxCut,
#              cuts.cBoxCut,
#              cuts.cBoxCut,
#              cuts.cNRecoCut,
#              cuts.cNRecoCut]

# #Cut keys
# cut_keys = ['Etheta',
#             'ele.len',
#             'ele.openangle',
#             'ele.electron',
#             'ele.dedx',
#             'ele.Etheta',
#             'ltrk.len',
#             'ltrk.len']

# #Bool keys - says whether event passes or fails cut
# bool_keys = ['cEtheta2',
#              'celeLen',
#              'celeOpenAngle',
#              'celeElectron',
#              'celededx',
#              'celeEtheta2',
#              'cTrkLenNShwE1',
#              'cTrkLenNShwG1']

# #Equalities - use greater than, less than, or between. For between cut ranges needs two lists
# equalities = ['lessthan',
#               'lessthan',
#               'lessthan',
#               'greaterthan',
#               'lessthan',
#               'lessthan',
#               'lessthan',
#               'lessthan']

# #nreco_keys for Nreco cut.
# nreco_keys = [None,
#               None,
#               None,
#               None,
#               None,
#               None,
#               'nshw',
#               'nshw']

# #nreco values for NReco cut
# nreco = [None,
#               None,
#               None,
#               None,
#               None,
#               None,
#               1,
#               1]

# #nreco_equalites values for NReco cut
# nreco_equalities = [None,
#               None,
#               None,
#               None,
#               None,
#               None,
#               'equal',
#               'greaterthan']

cut_order = [7,6,4,2,1,3,5] #Order for cuts to be applied
#cut_order = [4,2,1,3,5] #Order for cuts to be applied
#logging.info(f'{len(cut_ranges)},{len(cut_funcs)},{len(cut_keys)},{len(bool_keys)},{len(equalities)},{len(nreco_keys)},{len(nreco)},{len(nreco_equalities)}')

#Get no cut values
intsig_events = np.count_nonzero(events.evt_type == r'$\nu + e$') #Initial signal events
tot_events = len(events) #Initial signal events
pur,eff,f1,sig_events,tot_events = cuts.efficiency_purity_f1_calc(events,np.full(len(events),True),intsig_events=intsig_events)
sig_back_events = [intsig_events,tot_events] #signal,background events
modified_f1_min = [0,1e10] #cut_val,f1_val
np.savetxt(f'{folder_name}/effs/no_cuts.txt',[[0,eff],[0,pur],modified_f1_min,sig_back_events])
print(f'No cut, sig = {sig_back_events[0]}, back = {sig_back_events[1]}, cut val = None')

"""
Selection dictionary has the cut key as the key and in it is
cut_range : the range of cut values to test
cut_function : what type of cut
cut_key : which key it is in the df
equality : cut on values lessthan, greaterthan, equal?

Next 3 are for NReco cut only, set to None if using a box cut
nreco_key : which key to check for nreco objects, i.e. nshw
nreco : number of reco objects, i.e. nshw = 1 would be 1
nreco_equality : check for reco objects lessthan, greaterthan, equal, i.e. nshw = 1 would be equal
"""

for bool_key,sel_dict in nue_sel_dict.items():
  print(bool_key)
  #Extract values
  cut_vals = sel_dict['cut_range']
  cut_func = sel_dict['cut_function']
  cut_key = sel_dict['cut_key']
  equality = sel_dict['equality']
  nreco_key = sel_dict['nreco_key']
  nreco = sel_dict['nreco']
  nreco_equality = sel_dict['nreco_equality']
  #Efficiency,purity and f1 initialization
  effs = np.zeros(cut_vals.shape[0])
  purs = effs.copy()
  f1s = effs.copy()
  modified_f1s = effs.copy() #Modify so that max pur and eff are included in calc
  masks = np.full((cut_vals.shape[0],len(events)),True) #Temporarily store all masks, save best one based on modified f1 score

  #Find max eff,pur minimum f1
  f1_min = [0,1e10] #cut_val,f1_val
  modified_f1_min = [0,1e10] #cut_val,f1_val
  eff_max = [0,0] #cut_val,eff_val
  pur_max = [0,0] #cut_val,pur_val
  
  for i,cut_val in enumerate(cut_vals):
    #Apply cut from list of cut functions
    if cut_func == cuts.cBoxCut:
      events = cut_func(events,
                        cut_key=cut_key,
                        bool_key=bool_key,
                        cut_val=cut_val,
                        equality=equality)
      xlabel = f'{cut_key} cut value (keep events {equality})'
      save_name = f'eff_pur_f1_{bool_key}'
      title = None
    elif cut_func == cuts.cNRecoCut:
      events = cut_func(events,
                        cut_key=cut_key,
                        nreco_key=nreco_key,
                        nreco=nreco,
                        cut_val=cut_val,
                        bool_key=bool_key,
                        equality=equality,
                        nreco_equality=nreco_equality)
      xlabel = f'{cut_key} cut value (keep events {equality})'
      if nreco_equality == 'equal':
        title = f'{nreco_key} = {nreco}'
      elif nreco_equality == 'greaterthan':
        title = f'{nreco_key} > {nreco}'
      elif nreco_equality == 'lessthan':
        title = f'{nreco_key} > {nreco}'
      save_name = f'eff_pur_f1_{bool_key}_{nreco_key}_{nreco_equality}_{nreco}'
    #print(events.head())
    mask = events.loc[:,bool_key].values #get events that pass cut
    masks[i] = mask
    pur,eff,f1,sig_events,tot_events = cuts.efficiency_purity_f1_calc(events,mask,intsig_events=intsig_events)
    effs[i] = eff
    purs[i] = pur
    f1s[i] = f1

    #Assign values to max/min
    if eff > eff_max[1]:
      eff_max = [cut_val,eff]
    if pur > pur_max[1]:
      pur_max = [cut_val,pur]
    if f1 < f1_min[1] and f1 != 0: #second condition for dummy values
      f1_min = [cut_val,f1]
  modified_f1s = cuts.f1_calc(purs,effs,max_pur=max(purs),max_eff=max(effs)) #Get max pur+eff weighted f1 scores
  #logging.info(f'{modified_f1s}')
  #logging.info(f'{effs}')
  for k,f1 in enumerate(modified_f1s):
    if f1 < modified_f1_min[1] and f1 != 0:
      modified_f1_min = [cut_vals[k],f1]
      best_mask = masks[k] #Marks which events to keep and cut
      sig_back_events = [np.count_nonzero(events[best_mask].evt_type == r'$\nu + e$'),
                         np.count_nonzero(events[best_mask].evt_type != r'$\nu + e$')] #signal,background events
  
  np.savetxt(f'{folder_name}/effs/{save_name}.txt',[eff_max,pur_max,modified_f1_min,sig_back_events])
  if apply_cuts:
      events = events[best_mask] #Apply cuts
      tot_events = len(events)
  print(f'-- sig = {sig_back_events[0]:.2e}, back = {sig_back_events[1]:.2e}, cut val = {modified_f1_min[0]}')
  if make_plots:
    fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)
    ax2 = ax.twinx()
    ln1 = ax.plot(cut_vals,effs,label='Eff',color='blue')
    ax.scatter(eff_max[0],eff_max[1],color='blue')
    ln2 = ax.plot(cut_vals,purs,label='Pur',color='red')
    ax.scatter(pur_max[0],pur_max[1],color='red')
    ln3 = ax2.plot(cut_vals,modified_f1s,label='F1',color='green')
    ax2.scatter(modified_f1_min[0],modified_f1_min[1],color='green')

    ax.set_xlabel(xlabel)
    ax.set_ylabel('Pur and Eff')
    ax2.set_ylabel(r'F1 = $\frac{1}{2} ($max(pur)$/$pur$ + $max(eff)$/ $eff$ )$')
    ax.set_title(title)

    plotters.set_style(ax)
    plotters.set_style(ax2)

    #combine legends
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    lines = lines1 + lines2
    labels = labels1 + labels2
    ax.legend(lines, labels)
    plotters.save_plot(save_name,folder_name=folder_name+'/effs')
    ax.set_yscale('log')
    ax2.set_yscale('log')

    #plotters.save_plot(f'{save_name}_logy',folder_name=folder_name+'/effs')








