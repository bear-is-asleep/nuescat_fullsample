import pandas as pd
import seaborn as sns
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import logging

import sys
sys.path.append('../')
sys.path.append('/sbnd/app/users/brindenc/mypython/bc_utils')
from utils import plotters,pic
from CAFdata import *
import helpers
from datetime import date

#logging.basicConfig(level=logging.DEBUG)

day = date.today().strftime("%Y_%m_%d")
plotters.use_science_style()
#NuEScat_dir = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_6_fullsample_few_recocut_recoslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_8_fullsample_nue_few_recocut_recoslc'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_2_23_pure_nue_truth_cuts'
#NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_14_fullsample_softetheta_recoslc'
NuEScat_dir = '/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_29_fullsample_nue_trueav'
surname = ''
fnames = [f'cut_events_basic{surname}.root',f'cut_events_trk{surname}.root',f'cut_events_shw{surname}.root']
folder_name=f'{NuEScat_dir}/plots/{day}'
SAMPLE_POT = 3.08415e+20
NOM_POT = 10e20
NORMALIZATION = NOM_POT/SAMPLE_POT #Renormalize all counts to this number

colors = helpers.colors
categories = helpers.categories
def key_mask_values(x,key):
  """
  Mask extreme values such that we view important regions of the plots.
  These should absolutely not cut any events

  I realize now that we should not place cuts, but restrict x-axis (and save as seperate plot)
  """
  if key in ['true_Etheta','Etheta']:
    x = x[x<0.1]
  if key[-3:] == 'eng':
    x = x[x<4]
  if key[-4:] == 'dedx':
    x = x[x<10]
  if key[-3:] == 'gap':
    x = x[x<450]
  if key[-4:] == 'dens':
    x = x[x<40]
  if key[-3:] == 'len':
    x = x[x<300]
  return x

def key_set_lim(ax,key,axis='x'):
  """
  I realize now that we should not place cuts, but restrict x-axis (and save as seperate plot)
  """
  if axis == 'x':
    #if key in ['true_Etheta','Etheta']:
    #  ax.set_xlim([-0.01,0.1])
    if key[-3:] == 'eng':
      ax.set_xlim([-0.1,3])
    if key[-4:] == 'dedx':
      ax.set_xlim([-0.5,10])
    if key[-3:] == 'gap':
      ax.set_xlim([None,50])
    if key[-4:] == 'dens':
      ax.set_xlim([-5,40])
    if key[-3:] == 'len':
      ax.set_xlim([-10,300])
    if key[-9:] == 'openangle':
      ax.set_xlim([0,90])
  elif axis=='y':
    #if key in ['true_Etheta','Etheta']:
    #  ax.set_ylim([-0.01,0.1])
    #if key[-3:] == 'eng':
    #  ax.set_ylim([-0.1,4])
    if key[-4:] == 'dedx':
      ax.set_ylim([-0.5,10])
    if key[-3:] == 'gap':
      ax.set_ylim([None,50])
    if key[-4:] == 'dens':
      ax.set_ylim([-5,40])
    if key[-3:] == 'len':
      ax.set_ylim([-10,300])
    if key[-9:] == 'openangle':
      ax.set_ylim([0,90])
  return ax

def plot_df_keys(df,ftype='.png',keys_x=None,keys_y=None,plot_scat=False,plot_kde=True,
                 mask=False,**kwargs):
  """
  plot_df_keys(df)

  Parameters:
  df (pandas dataframe): dataframe that contains the data to be plotted
  keys_x: keys to plot on x axis
  keys_y: keys to plot on y axis
  plot_scat: True to plot scatters
  plot_kde: True to plot kde
  mask: apply a mask?

  Returns:
  None
  """
  counter = 0
  fnames = []
  alphas = [1,0.4,0.4,0.4,0.4,0.4]
  if keys_x is None: #Filter to just keys you want if necessary
    keys_x = df.keys()
  if keys_y is None: #Filter to just keys you want if necessary
    keys_y = df.keys()
  os.system(f'mkdir -p {folder_name}')
  for key_x in keys_x:
    for key_y in keys_y:
      #Conditions to avoid making too many plots
      if (key_x != key_y or key_x == 'evt_type') and not plot_scat: continue #skip scatter plots
      if (key_x == key_y) and not plot_kde: continue #skip kde plots
      if key_y == 'evt_type': continue #dont do evt type on y axis
      if key_x[0] == 'c' or key_y[0] == 'c':continue #skip cut
      xs = [] #x axis points
      ys = [] #y axis points
      labels = []
      alphas = np.zeros(len(event_names))
      fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)
      for i,name in enumerate(event_names):
        if name == r'$\nu + e$':
          alphas[i] = 0.9
        else: 
          if key_x == key_y:
            alphas[i] = 0.3
          else:
            alphas[i] = 0.1
        x = df.loc[df['evt_type']==name,key_x].values #get values
        if mask: #Apply masks for certain values
          x = key_mask_values(x,key_x)
        
        x,drop_ind = helpers.remove_dummy_values(x,
                                        dummy_val_list=[-9999,-999,-5,9999,999],
                                        is_df=False,
                                        return_drop_indeces=True)
        
        y = np.delete(df.loc[df['evt_type']==name,key_y].values,drop_ind) #use indeces from x for same size scatter plot
        xs.append(x)
        ys.append(y)
        labels.append(f'{name} ({int(NORMALIZATION* len(xs[i])):,})')
        if key_x != key_y:
          sns.scatterplot(x=xs[i],y=ys[i],label=labels[i],ax=ax,alpha=alphas[i],**kwargs)
          ax.set_xlabel(f'{key_x}')
          ax.set_ylabel(f'{key_y}')
          prefix = 'scat'
        if key_x == key_y and key_x != 'evt_type':
          #Set bin widths, fine tuning
          # if key_x in ['Etheta','true_Etheta'] and name != r'$\nu + e':
          #   bw = 0.01
          #   ax = sns.kdeplot(x=xs[i],ax=ax,label=labels[i],warn_singular=False,alpha=alphas[i],bw_method=bw,**kwargs)
          if key_x in ['nreco','nshw','ntrk','nslc','nstub']:
            bw = 0.5
            ax = sns.kdeplot(x=xs[i],ax=ax,label=labels[i],warn_singular=False,alpha=alphas[i],bw_method=bw,**kwargs)
          else:
            ax = sns.kdeplot(x=xs[i],ax=ax,label=labels[i],warn_singular=False,alpha=alphas[i],**kwargs)
          ax.set_xlabel(f'{key_x}')
          ax.set_ylabel('Counts')
          prefix = 'dist'
      ax.legend()
      if key_x[-6:] != 'Etheta':
        key_set_lim(ax,key_x)
      plotters.set_style(ax,legend_size=12)
      fnames.append(f'{prefix}_{key_x.replace(".","_")}_{key_y.replace(".","_")}')
      plt.savefig(f'{fnames[counter]}{ftype}',bbox_inches = "tight",dpi=300)
      print(fnames[counter])
      plt.clf()
      plt.close('all')
      del ax
      del fig
      counter+=1
  for fname in fnames:
    os.system("mv " + fname + ftype + f' {folder_name}/')

def hist2d_true_reco(reco,true,reco_label,true_label,title=None,plot_line=False,label_boxes=False,
                   **kwargs):
  """
  Plot 2d hist with x axis as reco, y axis as true
  plot_line to plot y=x line
  """
  lower = np.min([0,min(reco),min(true)])
  upper = np.max([max(reco),max(true)])

  fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)

  if plot_line:
    xy = [lower,upper]
    ax.plot(xy,xy,ls='--',color='red')
    ax.set_xlim([lower,upper])
    ax.set_ylim([lower,upper])
    hist,xbins,ybins,im = ax.hist2d(reco,true,range=[xy,xy],**kwargs)
  else:
    hist,xbins,ybins,im = ax.hist2d(reco,true,**kwargs)
  if label_boxes:
    for i in range(len(ybins)-1):
      for j in range(len(xbins)-1):
        ax.text(xbins[j]+0.5,ybins[i]+0.5, f'{hist.T[i,j]:.0f}', 
                color="w", ha="center", va="center", fontweight="bold",fontsize=16)
  ax.set_xlabel(f'{reco_label}')
  ax.set_ylabel(f'{true_label}')
  if title is not None:
    ax.set_title(title)
  plotters.set_style(ax)
  return fig,ax

def plot_true_reco_err(reco,true,reco_label,true_label,title=None,**kwargs):
  """
  plot histogram with error between true and reco
  """
  err = (reco-true)/true
  fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)

  ax.hist(err,**kwargs)
  ax.set_xlabel(f'({reco_label}-{true_label})/{true_label}')
  if title is not None:
    ax.set_title(title)
  plotters.set_style(ax)
  return fig,ax

def plot_hist(df,key_x,title=None,**kwargs):
  """
  plot a histogram with several bins
  """
  labels = []
  xs = []
  weights = []
  for i,name in enumerate(event_names):
    x = df.loc[df['evt_type']==name,key_x].values
    #Mask certain values
    x = helpers.remove_dummy_values(x,
                                    dummy_val_list=[-9999,-999,-5,9999,999,0,-1],
                                    is_df=False)
    if len(x) == 0: continue #Skip event types that have no events
    xs.append(x)
    weights.append(np.full(x.shape,NORMALIZATION))
    labels.append(f'{name} ({int(NORMALIZATION* len(x)):,})')
  fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)
  ax.hist(xs,label=labels,edgecolor=colors[i],stacked=True,weights=weights,**kwargs)
  ax.set_xlabel(key_x)
  ax.set_ylabel('Count')
  if title is not None:
    ax.set_title(title)
  ax.legend()
  plotters.set_style(ax)
  return fig,ax

def plot_masked_vals(df,x_key,masks,labels,use_kde=True,y_key=None,title=None,**kwargs):
  """
  df contains dataframe with keys
  Plots scatter and histograms using masked values and labels them.
  x is the data array in x[0] and is the label in x[1]
  y is optional and will make this a scatter plot

  """
  xs = []
  cnt_labels = []
  colors = plotters.get_colors_from_map(len(labels),cmap='nipy_spectral') #Get sequential colors
  x = df.loc[:,x_key].values
  logging.debug(f'{len(x)},{len(masks[0])}')
  x,drop_inds = helpers.remove_dummy_values(x,
                                  dummy_val_list=[-9999,-999,-5,9999,999],
                                  is_df=False,
                                  return_drop_indeces=True)
  if y_key is not None:
    y = np.delete(df.loc[:,y_key].values,drop_inds) #use indeces from x for same size scatter plot
    y,drop_indsy = helpers.remove_dummy_values(y,
                                  dummy_val_list=[-9999,-999,-5,9999,999],
                                  is_df=False,
                                  return_drop_indeces=True)
    x = np.delete(x,drop_indsy)
    drop_inds = np.concatenate([drop_inds,drop_indsy])
  masks = [np.delete(mask,drop_inds) for mask in masks] #Drop inds that are dummy values
  fig,ax = plt.subplots(figsize=(8,6),tight_layout=True)
  for i,mask in enumerate(masks):
    logging.debug(f'{len(x),len(mask)},{len(drop_inds)},{len(np.unique(drop_inds))}')
    #mask = np.delete(mask,drop_inds)
    logging.debug(f'{len(x),len(mask)}')
    #if len(x) != len(mask): continue #skip buggy events, not sure why these happen
    if len(x[mask]) == 0: continue #skip empty events
    label = f'{labels[i]} ({len(x[mask])})'
    xs.append(x[mask])
    cnt_labels.append(label)
    if y_key is not None:
      ax.scatter(x[mask],y[mask],label=label,color=colors[i],**kwargs)
    elif use_kde:
      ax = sns.kdeplot(x=x[mask],ax=ax,label=label,color=colors[i],warn_singular=False,**kwargs)
  else:
    if len(xs) == 0:
      ax.hist(xs,label=cnt_labels,**kwargs)
    else:  
      ax.hist(xs,label=cnt_labels,color=colors[:len(xs)],**kwargs)
  if title is not None:
    ax.set_title(title)
  if y_key is not None:
    ax.set_ylabel(y_key)
  else:
    ax.set_ylabel('Prob.')
  ax.set_xlabel(x_key)
  ax.legend()
  plotters.set_style(ax,legend_size=12)
  #bbox_to_anchor=(1.05,1)
  return fig,ax

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

#events = pd.read_csv(f'{NuEScat_dir}/cut_events.csv')
#events = events.set_index(['run','subrun','evt'])
#print(events.loc[:,('ntrk','nshw')])
print('-----Cleaning data')
#Mask event type and get their names
events = helpers.set_event_type(events)
events.loc[:,'nreco'] = events.nshw + events.ntrk
events.loc[:,'theta_err'] = (events.reco_theta-events.true_theta)/events.true_theta
events.loc[:,'Q2_nuecc'] = helpers.Q2_nuecc(helpers.remove_dummy_values(events.loc[:,'true_slice_eng'].values,dummy_val_list=[-9999,-999,999,9999,-1,0])
                                            ,helpers.remove_dummy_values(events.loc[:,'true_theta']))
#events = cuts.apply_cuts(events)
event_names = categories

print('-----Seperate data types')
#Isolate data from different reconstructed objects, shws trks
lshw_keys = [key for key in events.keys() if key[:4] == 'lshw']
lshw_keys.extend(['evt_type'])
slshw_keys = [key for key in events.keys() if key[:5] == 'slshw']
slshw_keys.extend(['evt_type'])
ltrk_keys = [key for key in events.keys() if key[:4] == 'ltrk']
ltrk_keys.extend(['evt_type'])
sltrk_keys = [key for key in events.keys() if key[:5] == 'sltrk']
sltrk_keys.extend(['evt_type'])

#Get all events from dataframe
lshw = events.loc[:,lshw_keys]
slshw = events.loc[:,slshw_keys]
ltrk = events.loc[:,ltrk_keys]
sltrk = events.loc[:,sltrk_keys]

#Remove events filled with dummy values
# lshw = helpers.remove_dummy_values(lshw)
# slshw = helpers.remove_dummy_values(slshw)
# ltrk = helpers.remove_dummy_values(ltrk)
# sltrk = helpers.remove_dummy_values(sltrk)

other_keys = []
for key in events.keys():
  if key not in slshw_keys and key not in sltrk_keys and key not in ltrk_keys and key not in lshw_keys:
    if key != 'run' and key != 'subrun' and key != 'evt':
      other_keys.append(key)
other_keys.extend(['evt_type'])
other_events = events.loc[:,other_keys]
#events = helpers.remove_dummy_values(events,dummy_val_list=[-9999,-999,999,9999,-1,-5,0])

#Get reco and true labels and arrays
reco_labels_other = ['nshw','ntrk','reco_eng','reco_theta',
                     'reco_vtx.x','reco_vtx.y','reco_vtx.z','Etheta','nreco','nreco']
reco_array_other = [other_events.loc[:,lab] for lab in reco_labels_other]
true_labels_other = ['truenshw','truentrk','true_slice_eng','true_theta',
                     'true_vtx.x','true_vtx.y','true_vtx.z','true_Etheta','nshw','ntrk']
true_array_other = [other_events.loc[:,lab] for lab in true_labels_other]
#Append objects not in the array
true_array_other.append(other_events.true_slice_eng*other_events.true_theta**2)
true_array_other.append('true_Etheta2')

int_bins = np.arange(0,8,1)
eng_bins = np.arange(0,3,0.2)
theta_bins = np.arange(0,0.8,0.05)
x_bins = np.arange(-250,250,10)
z_bins = np.arange(-50,550,10)
etheta_bins = np.arange(0,stop=0.02,step=0.002)
err_bins = np.arange(-1,1,0.05)
theta_err_bins = np.arange(-2,26,1)
de_dx_bins = np.arange(0,10,0.5)
dens_bins = np.arange(0,50,2.5)
len_bins = np.arange(0,500,10)
bins_other = [int_bins,int_bins,eng_bins,theta_bins,x_bins,x_bins,z_bins,
              etheta_bins,int_bins,int_bins]

#Masks
# ltrk_electron_mask = events.loc[:,'ltrk.true.pdg'].values == 11
# sltrk_electron_mask = events.loc[:,'sltrk.true.pdg'].values == 11
# lshw_electron_mask = events.loc[:,'lshw.true.pdg'].values == 11
# slshw_electron_mask = events.loc[:,'slshw.true.pdg'].values == 11

# ltrk_muon_mask = abs(events.loc[:,'ltrk.true.pdg'].values) == 13
# sltrk_muon_mask = abs(events.loc[:,'sltrk.true.pdg'].values) == 13
# lshw_muon_mask = abs(events.loc[:,'lshw.true.pdg'].values) == 13
# slshw_muon_mask = abs(events.loc[:,'slshw.true.pdg'].values) == 13

# masks = [ltrk_electron_mask,sltrk_electron_mask,lshw_electron_mask,slshw_electron_mask,
#          ltrk_muon_mask,sltrk_muon_mask,lshw_muon_mask,slshw_muon_mask]
# labels = ['ele ltrk','ele sltrk','ele lshw','ele slshw',
#           'muon ltrk','muon sltrk','muon lshw','muon slshw']
# masks_ele = masks[:4]
# labels_ele = labels[:4]
# masks_muon = masks[4:]
# labels_muon = labels[4:]

#nreco masks
ntrk0_mask = events.loc[:,'ntrk'].values == 0
ntrk1_mask = events.loc[:,'ntrk'].values >= 1
#ntrk2_mask = events.loc[:,'ntrk'].values >= 2

nshw0_mask = events.loc[:,'nshw'].values == 0
nshw1_mask = events.loc[:,'nshw'].values == 1
nshw2_mask = events.loc[:,'nshw'].values >= 2

#Combine masks
ntrk0_nshw0_mask = np.logical_and(ntrk0_mask,nshw0_mask)
ntrk1_nshw0_mask = np.logical_and(ntrk1_mask,nshw0_mask)
ntrk0_nshw1_mask = np.logical_and(ntrk0_mask,nshw1_mask)
ntrk1_nshw1_mask = np.logical_and(ntrk1_mask,nshw1_mask)
ntrk0_nshw2_mask = np.logical_and(ntrk0_mask,nshw2_mask)
ntrk1_nshw2_mask = np.logical_and(ntrk1_mask,nshw2_mask)

masks = [ntrk0_nshw0_mask,
        ntrk1_nshw0_mask,
        ntrk0_nshw1_mask,
        ntrk1_nshw1_mask,
        ntrk0_nshw2_mask,
        ntrk1_nshw2_mask]
mask_labels = ['ntrk = 0 & nshw = 0',
               'ntrk > 0 & nshw = 0',
               'ntrk = 0 & nshw = 1',
               'ntrk > 0 & nshw = 1',
               'ntrk = 0 & nshw > 1',
               'ntrk > 0 & nshw > 1']

plot_mask_keys = ['reco_eng','true_theta','reco_theta','Etheta','true_Etheta','true_slice_eng','nstub','theta_err']

print('-----Test plots')
# x_key = 'lshw.electron'
# fig,ax = plot_masked_vals(events,x_key,masks,mask_labels)
#plotters.save_plot(f'pdgmask_ele_{x_key}',folder_name=folder_name,fig=fig)
#plot_df_keys(other_events,keys_x=['reco_eng'])
# plot_df_keys(lshw,keys_x=['lshw.len','lshw.cnvgap','lshw.dens'])
if False:
  for rl,tl,ra,ta,bins in zip(reco_labels_other,true_labels_other,reco_array_other,true_array_other,bins_other):
      if rl == 'reco_theta':
        fig,ax = plot_true_reco_err(ra,ta,rl,tl,title=r'$\nu+e$',bins=theta_err_bins,alpha=0.3,linewidth=3)
      else:
        fig,ax = plot_true_reco_err(ra,ta,rl,tl,title=r'$\nu+e$',bins=err_bins,alpha=0.3,linewidth=3)
      ax.axvline(0,ls='--',color='red')
      plotters.save_plot(f'err_{tl}_{rl}',folder_name=folder_name,fig=fig)
      plt.close('all')

# fig,ax = plot_hist(events,'reco_eng',alpha=0.5,linewidth=2,bins=eng_bins)
# plotters.save_plot(f'hist_Etheta_zoom',folder_name=folder_name,fig=fig)

# plot_df_keys(lshw,keys_x=['lshw.eng'],plot_kde=False,plot_scat=True)

print('-----Starting plots')
#Truth plots
fig,ax = plot_hist(events,'true_Etheta',bins=np.arange(0,1,0.025))
ax.set_xlabel(r'$E_e\theta_e^2$ [Gev rad$^2$]')
plotters.save_plot('aps_etheta2_singleele',folder_name=folder_name)

fig,ax = plot_hist(events,'true_Etheta',bins=np.arange(0,0.025,0.001))
ax.set_xlabel(r'$E_e\theta_e^2$ [Gev rad$^2$]')
#ax.axvline(x=2*0.511*1e-3,linestyle='--')
#ax.text(2*0.511*1e-3+0.001,200,r'$2m_e$')
plotters.save_plot('aps_etheta2_singleele_zoom',folder_name=folder_name)

fig,ax = plot_hist(events,'true_slice_eng',bins=np.arange(0,3,0.1))
ax.set_xlabel(r'$E_e$ [Gev]')
plotters.save_plot('aps_sliceeng_singleele',folder_name=folder_name)
ax.set_yscale('log')
plotters.save_plot('aps_sliceeng_singleele_logy',folder_name=folder_name)

fig,ax = plot_hist(events,'true_theta',bins=np.arange(0,1,0.05))
ax.set_xlabel(r'$\theta_e$ [rad]')
plotters.save_plot('aps_theta_singleele',folder_name=folder_name)
ax.set_yscale('log')
ax.set_ylim([None,4e3])
plotters.save_plot('aps_theta_singleele_logy',folder_name=folder_name)

fig,ax = plot_hist(events,'Q2_nuecc',bins=np.arange(0,1,0.05))
ax.set_xlabel(r'$Q^2$ [GeV$^2$] ($\nu_e$ CC)')
plotters.save_plot('aps_Q2_singleele',folder_name=folder_name)
ax.set_yscale('log')
ax.set_ylim([None,4e3])
plotters.save_plot('aps_Q2_singleele_logy',folder_name=folder_name)

plt.close('all')

#Make comps to true type of second 
# for i,x_key in enumerate(events):
#   for j,y_key in enumerate(events):
#     print(f'mask plot {x_key} {y_key}')
#     if x_key == y_key:
#       fig,ax = plot_masked_vals(events,x_key,masks_ele,labels_ele,alpha=0.8)
#       key_set_lim(ax,x_key)
#       pname = f'dist_ele__{x_key.replace(".","_")}'
#     else:
#       fig,ax = plot_masked_vals(events,x_key,masks_ele,labels_ele,use_kde=False,y_key=y_key,
#                                 alpha=0.5)
#       key_set_lim(ax,x_key)
#       key_set_lim(ax,y_key,axis='y')
#       pname = f'scat_ele__{x_key.replace(".","_")}_{y_key.replace(".","_")}'
#     plotters.save_plot(pname,folder_name=folder_name+'/pdgmask',fig=fig)
#     plt.close('all')

#Make comps to true type
if False:
  for i,x_key in enumerate(events):
    for j,y_key in enumerate(events):
      #if x_key not in plot_mask_keys: continue
      if y_key in ['evt_type']: continue
      
      if x_key == y_key and x_key not in ['evt_type']:
        print(f'mask plot {x_key} {y_key}')
        #fig,ax = plot_masked_vals(events,x_key,masks,mask_labels,alpha=0.7)
        #key_set_lim(ax,x_key)
        #pname = f'dist_{x_key.replace(".","_")}'
      else:
        continue
        fig,ax = plot_masked_vals(events,x_key,masks,mask_labels,use_kde=False,y_key=y_key,
                                  alpha=0.3)
        #key_set_lim(ax,x_key)
        #key_set_lim(ax,y_key,axis='y')
        pname = f'scat_{x_key.replace(".","_")}_{y_key.replace(".","_")}'
      #plotters.save_plot(pname,folder_name=folder_name+'/pdgmask',fig=fig)
      #ax.cla()
      #plt.clf()
      #plt.close('all')
      if x_key == y_key and x_key not in ['evt_type']:
        if x_key in ['truentrk','truenshw','nele']: continue
        elif x_key in ['nshw','nreco','ntrk','nstub','truenshw','truentrk','nele']:
          bins = np.arange(0,max(events.loc[:,x_key]),1) #arange bins into integer numbers
        elif x_key[-3:] == 'eng':
          bins=eng_bins
        elif x_key[-4:] in ['vgap','dens']:
          bins=dens_bins
        elif x_key[-3:] in ['len']:
          bins=len_bins
        elif x_key[-6:] in ['Etheta']:
          bins=etheta_bins
        elif x_key[-5:] in ['theta']:
          bins = theta_bins
        elif x_key in ['theta_err']:
          bins = theta_err_bins
        elif x_key[-4:] in ['dedx']:
          bins = de_dx_bins
        else: bins = 20
        fig,ax = plot_masked_vals(events,x_key,masks,mask_labels,use_kde=False,
                                  alpha=0.7,bins=bins,histtype='step',linewidth=2)
        key_set_lim(ax,x_key)
        # if x_key[-3:] == 'eng':
        #   ax.set_ylim([0,180])
        # if x_key[-5:] == 'theta':
        #   ax.set_ylim([0,310])
        if x_key == 'lshw.photon':
          ax.set_ylim([0,55])
        pname = f'hist_{x_key.replace(".","_")}'
        plotters.save_plot(pname,folder_name=folder_name+'/pdgmask',fig=fig)
        ax.set_yscale('log')
        plotters.save_plot(pname+'_logy',folder_name=folder_name+'/pdgmask',fig=fig)
        plt.close('all')


#Make histograms
if False:
  for i,key in enumerate(events):
    if key == 'evt_type': continue

    # if key[-4:] in ['dedx']:
    #   fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=de_dx_bins)
    if key in ['nshw','nreco','ntrk','nstub','truenshw','truentrk','nele']:
      bins = np.arange(0,max(events.loc[:,key]),step=1) #arange bins into integer numbers
    #   fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=bins)
    # elif key[-3:] == 'eng':
    #   fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=eng_bins)
    # elif key[-4:] in ['vgap','dens']:
    #   fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=dens_bins)
    # elif key[-3:] in ['len']:
    #   fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=len_bins)
    else: bins = 10
    fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=20)
    # if key != 'Etheta':
    #   key_set_lim(ax,key)
    key = key.replace('.','_')
    plotters.save_plot(f'hist_{key}',folder_name=folder_name,fig=fig)
    # if key[-6:] == 'Etheta':
    #   fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=etheta_bins)
    #   key = f'{key}_zoom'
    #   plotters.save_plot(f'hist_{key}',folder_name=folder_name,fig=fig)
    #ax.set_yscale('log')
    #plotters.save_plot(f'hist_{key}_logy',folder_name=folder_name,fig=fig)
    plt.close('all')
    print(key)

#plot_df_keys(other_events)
#plot_df_keys(lshw)
#plot_df_keys(slshw)
#plot_df_keys(ltrk)
#plot_df_keys(sltrk)
scat_keys_other = ['reco_eng','reco_theta','nshw','ntrk','nstub']
#plot_df_keys(other_events,keys_x=scat_keys_other,plot_kde=False,plot_scat=True)


#Make plots
if False:
  for rl,tl,ra,ta,bins in zip(reco_labels_other,true_labels_other,reco_array_other,true_array_other,bins_other):
    if rl in ['nshw','ntrk','nreco']:
      fig,ax = hist2d_true_reco(ra,ta,rl,tl,title=r'$\nu+e$',bins=bins,cmap='viridis_r',label_boxes=True)
      plotters.save_plot(f'hist2d_{tl}_{rl}',folder_name=folder_name,fig=fig)
    else:
      fig,ax = hist2d_true_reco(ra,ta,rl,tl,title=r'$\nu+e$',plot_line=True,bins=bins)
      plotters.save_plot(f'hist2d_{tl}_{rl}',folder_name=folder_name,fig=fig)
      if rl == 'reco_theta':
        fig,ax = plot_true_reco_err(ra,ta,rl,tl,title=r'$\nu+e$',bins=theta_err_bins,alpha=0.3,linewidth=3)
      else:
        fig,ax = plot_true_reco_err(ra,ta,rl,tl,title=r'$\nu+e$',bins=err_bins,alpha=0.3,linewidth=3)
      ax.axvline(0,ls='--',color='red')
      plotters.save_plot(f'err_{tl}_{rl}',folder_name=folder_name,fig=fig)
    plt.close('all')
