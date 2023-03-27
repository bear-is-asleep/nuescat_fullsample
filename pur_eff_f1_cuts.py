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
make_cut_plots = True
apply_cuts = False
save_final_selection = False
stride = 1 #select every "stride" events
nue_dict = nue_search_dict #Select if you want to search for or use parameters


day = date.today().strftime("%Y_%m_%d")
print(day)
plotters.use_science_style()
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_2_27_fullsample_nue_no_cuts_truthslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_1_fullsample_few_recocut_recoslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_9_fullsample_few_recocut_recoslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_13_fullsample_nue_few_recocut_recoslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_14_pure_nue_truefv_recoslc'
NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_14_fullsample_softetheta_recoslc'
folder_name = f'{NuEScat_dir}/cut_search_{day}'
info_name = f'{folder_name}/effs' #Save txt files
os.system(f'mkdir -p {folder_name}')
os.system(f'mkdir -p {info_name}')

surname = ''
fnames = [f'cut_events_basic{surname}.root',f'cut_events_trk{surname}.root',f'cut_events_shw{surname}.root']
#folder_name=f'{NuEScat_dir}/plots'

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
events = events.sample(n=int(len(events)/stride))

print('-----Cleaning data')
#Mask event type and get their names
events = helpers.set_event_type(events)
events.loc[:,'nreco'] = events.nshw + events.ntrk
#Calc etheta
events = helpers.calc_etheta_reco(events,'ltrk.')
events = helpers.calc_etheta_reco(events,'sltrk.')
events = helpers.calc_etheta_reco(events,'lshw.')
events = helpers.calc_etheta_reco(events,'slshw.')

#Options to mask out cosmics, or select electron based on clustering
events = cuts.mask_cosmics(events,time_threshold=[0.,1.6]) #Mask all events that have crt activity
#events = cuts.select_electron(events) #Should just be lshw

print('-----Make efficiency plots')


#Get no cut values
intsig_events = np.count_nonzero(events.evt_type == r'$\nu + e$') #Initial signal events
back_events = np.count_nonzero(events.evt_type != r'$\nu + e$') #Initial background events
pur,eff,f1,sig_events,tot_events = cuts.efficiency_purity_f1_calc(events,np.full(len(events),True),intsig_events=intsig_events)
sig_back_events = [intsig_events,back_events] #signal,background events
modified_f1_min = [0,1e10] #cut_val,f1_val
np.savetxt(f'{info_name}/no_cuts.txt',[[0,eff],[0,pur],modified_f1_min,sig_back_events])
print(f'No cut, sig = {sig_back_events[0]:.2e}, back = {sig_back_events[1]:.2e}, cut val = None')

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

#Store event count of each type after each cut
event_counts = np.zeros((len(nue_dict)+1,len(helpers.categories)))
cut_labels = ['']*(len(nue_dict)+1) #Store cut names in list
eff_overall = np.ones(len(nue_dict)+1) #Save efficiency values
pur_overall = np.zeros(len(nue_dict)+1) #Save all purity values
cut_count = 0

#Initialize no cut values
print('-- Remaining events : ')
for q,cat in enumerate(helpers.categories):
  cat_events = len(events[events.loc[:,'evt_type'] == cat])
  event_counts[cut_count][q] = cat_events
  print(f'---- {cat} = {cat_events:.2e}')
cut_labels[0] = 'No Cut'

#for bool_key,sel_dict in nue_dict.items():
for bool_key,sel_dict in nue_dict.items():
  cut_count += 1
  ratio_max = -1 #Find ratio between signal and background lost
  print(bool_key)

  #Extract dictionary values
  cut_vals = sel_dict['cut_range']
  if isinstance(cut_vals,list):
    cut_vals = np.array(cut_vals)
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
  ratios = effs.copy()
  ratios = effs.copy() #Modify so that max pur and eff are included in calc
  masks = np.full((cut_vals.shape[0],len(events)),True) #Temporarily store all masks, save best one based on modified f1 score

  #Find max eff,pur minimum f1
  f1_min = [0,1e10] #cut_val,f1_val
  modified_f1_min = [0,1e10] #cut_val,f1_val
  eff_max = [0,0] #cut_val,eff_val
  pur_max = [0,0] #cut_val,pur_val
  ratio_max_vals = [0,0] #cut_val,ratio_val
  
  #Iterate over all cut values, select best one based on ratio
  for i,cut_val in enumerate(cut_vals):
    #Apply cut from list of cut functions
    if cut_func == cuts.cBoxCut:
      events = cut_func(events,
                        cut_key=cut_key,
                        bool_key=bool_key,
                        cut_val=cut_val,
                        equality=equality,
                        use_abs=False)
      
      #Labels for pur,eff,f1 plots
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
      #Labels for pur,eff,f1 plots
      xlabel = f'{cut_key} cut value (keep events {equality})'
      if nreco_equality == 'equal':
        title = f'{nreco_key} = {nreco}'
      elif nreco_equality == 'greaterthan':
        title = f'{nreco_key} > {nreco}'
      elif nreco_equality == 'lessthan':
        title = f'{nreco_key} > {nreco}'
      save_name = f'eff_pur_f1_{bool_key}_{nreco_key}_{nreco_equality}_{nreco}'
    #get events that pass cut
    mask = events.loc[:,bool_key].values 
    masks[i] = mask

    #Extract purity,eff,f1 scores
    pur,eff,f1,sig_events,tot_events = cuts.efficiency_purity_f1_calc(events,mask,intsig_events=intsig_events)
    back_events = tot_events-sig_events
    effs[i] = eff
    purs[i] = pur
    f1s[i] = f1

    #Find percent of events lost for bckgnd and signal
    total_background_events = np.sum(event_counts[cut_count-1,1:]) #Previous cut background events
    total_signal_events = event_counts[cut_count-1,0] #Previous cut signal events
    if total_background_events == 0:
      lost_back = 0
    else:
      lost_back = (total_background_events-back_events)/total_background_events
    if total_signal_events == 0:
      lost_sig = 0
      ratio = 0
    else:
      lost_sig = (total_signal_events-sig_events)/total_signal_events
      if lost_sig == 0: 
        lost_sig = 2*total_signal_events #This is equivalent to losing 1/2 signal event, we should not punsish the ratio for not losing signal events
      #####This controls the selection ratio####
      ratio = lost_back**(1)/lost_sig**(1)
    #print('cut_val...',cut_val,ratio,lost_sig,sig_events,tot_events,lost_back,
    #      total_background_events,total_signal_events,event_counts)
      ##########################################
    ratios[i] = ratio
    #Assign values to max/min
    if eff > eff_max[1]:
      eff_max = [cut_val,eff]
    if pur > pur_max[1]:
      pur_max = [cut_val,pur]
    if f1 < f1_min[1] and f1 != 0: #second condition for dummy values
      f1_min = [cut_val,f1]
    if ratio > ratio_max: #ratio of background/signal loss defines best cut value
      ratio_max = ratio
      best_mask = masks[i] #Marks which events to keep and cut
      sig_back_events = [np.count_nonzero(events[best_mask].evt_type == r'$\nu + e$'),
                         np.count_nonzero(events[best_mask].evt_type != r'$\nu + e$')] #signal,background events
      eff_overall[cut_count] = eff
      pur_overall[cut_count] = pur
      if equality == 'lessthan':
        cut_labels[cut_count] = f'{cut_key} < {cut_val}'
      elif equality == 'greaterthan':
        cut_labels[cut_count] = f'{cut_key} > {cut_val}'
      elif equality == 'equal':
        cut_labels[cut_count] = f'{cut_key} = {cut_val}' #Save cut label
      ratio_max_vals = [cut_val,ratio_max]
  #Outside cut_val searches, calculate modified f1 score
  modified_f1s = cuts.f1_calc(purs,effs,max_pur=max(purs),max_eff=max(effs)) #Get max pur+eff weighted f1 scores
  for k,f1 in enumerate(modified_f1s):
    if f1 < modified_f1_min[1] and f1 != 0:
      modified_f1_min = [cut_vals[k],f1]

  #Save efficieny, purity, f1, and sig/back events for best cut
  np.savetxt(f'{info_name}/{save_name}.txt',[eff_max,pur_max,modified_f1_min,sig_back_events,ratio_max_vals])
  
  if apply_cuts:
      events = events[best_mask] #Apply cuts

  #Display signal and background events after cut_val
  print(f'-- sig = {sig_back_events[0]:.2e}, back = {sig_back_events[1]:.2e}, cut val = {modified_f1_min[0]}')
  print(f'-- pur = {pur_overall[cut_count]:.2f}, eff = {eff_overall[cut_count]:.2f}, ratio = {ratio_max:.2e}')
  if make_plots: #Make purity and eff. plots, best used for cut search
    fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)
    ax2 = ax.twinx()
    ln1 = ax.plot(cut_vals,effs,label='Eff',color='blue',ls='-.',alpha=0.7)
    ax.scatter(eff_max[0],eff_max[1],color='blue',alpha=0.7)
    ln2 = ax.plot(cut_vals,purs,label='Pur',color='red',alpha=0.7)
    ax.scatter(pur_max[0],pur_max[1],color='red',alpha=0.7)
    # ln3 = ax2.plot(cut_vals,modified_f1s,label='F1',color='green',ls='--',alpha=0.7)
    # ax2.scatter(modified_f1_min[0],modified_f1_min[1],color='green',alpha=0.7)
    ln3 = ax2.plot(cut_vals,ratios,label='F1',color='green',ls='--',alpha=0.7)
    ax2.scatter(ratio_max_vals[0],ratio_max_vals[1],color='green',alpha=0.7)

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
    plotters.save_plot(save_name,folder_name=folder_name)
    ax.set_yscale('log')
    ax2.set_yscale('log')
    plotters.save_plot(f'{save_name}_logy',folder_name=folder_name)
    plt.close('all')

  #Store information across cuts
  print('-- Remaining events : ')
  for q,cat in enumerate(helpers.categories):
    cat_events = len(events[events.loc[:,'evt_type'] == cat])
    event_counts[cut_count][q] = cat_events
    print(f'---- {cat} = {cat_events:.2e}')
  


if make_cut_plots: #Make cut selection plots
  fig,ax = plt.subplots(figsize=(10,6),tight_layout=True)
  ax2 = ax.twinx()
  for i,cat in enumerate(helpers.categories):
    if cat == r'$\nu + e$': continue #Plot on separate axis
    else:
      alpha = 0.4
    ax.plot(cut_labels,event_counts[:,i],label=cat,color=helpers.colors[i],ls='-',
            alpha=alpha)
  ax2.plot(cut_labels,event_counts[:,0],label=r'$\nu + e$',color='blue',alpha=0.9)
  ax.set_xticklabels(cut_labels, rotation=45, ha='right')
  ax.set_ylabel('Background Events')
  ax2.set_ylabel('Signal Events')
  ax.grid(True,axis='x')
  #combine legends
  lines1, labels1 = ax.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  lines = lines1 + lines2
  labels = labels1 + labels2
  plotters.set_style(ax,tick_size=12)
  plotters.set_style(ax2,tick_size=12)
  ax.legend(lines, labels)
  ax2.set_ylim([1,None])
  ax.set_ylim([1,None])

  plotters.save_plot('cuts',folder_name=folder_name)
  ax.set_yscale('log')
  ax2.set_yscale('log')
  plotters.save_plot('cuts_logy',folder_name=folder_name)

  #Pur,eff,% of signal/background cut
  fig,ax = plt.subplots(figsize=(8,6),tight_layout=True)
  ax2 = ax.twinx()
  ax.plot(cut_labels,eff_overall,label='eff',color='blue',alpha=0.7)
  ax.plot(cut_labels,pur_overall,label='pur',color='red',alpha=0.7)
  lost_signal = -100*np.diff(event_counts[:,0])/event_counts[1:,0]
  ax2.plot(cut_labels[1:],lost_signal,
           label=r'$\nu + e$',color='green',alpha=0.9)
  event_counts_back = np.sum(event_counts[:,1:],axis=1)
  lost_background = -100*np.diff(event_counts_back)/event_counts_back[1:]
  ax2.plot(cut_labels[1:],lost_background,
           label=r'Background',color='brown',alpha=0.6)
  ax2.plot(cut_labels[1:],lost_background/lost_signal,
           label=r'Ratio',color='cyan',alpha=0.6)
  ax.set_xticklabels(cut_labels, rotation=45, ha='right')
  ax.set_ylabel('Purity/Efficiency')
  ax2.set_ylabel('% of events lost')
  ax.grid(True,axis='x')
  #combine legends
  lines1, labels1 = ax.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  lines = lines1 + lines2
  labels = labels1 + labels2
  plotters.set_style(ax,tick_size=12)
  plotters.set_style(ax2,tick_size=12)
  ax.legend(lines, labels)
  plotters.save_plot('pur_eff_lostevents',folder_name=folder_name)
  ax.set_yscale('log')
  ax2.set_yscale('log')
  plotters.save_plot('pur_eff_lostevents_logy',folder_name=folder_name)
  plt.close('all')

if save_final_selection:
  events.to_csv(f'{NuEScat_dir}/cut_events.csv')
