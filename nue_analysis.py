import pandas as pd
import seaborn as sns
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import cuts

import sys
sys.path.append('../')
sys.path.append('/sbnd/app/users/brindenc/mypython/bc_utils')
from utils import plotters,pic
from CAFdata import *
import helpers
from datetime import date

day = date.today().strftime("%Y_%m_%d")
plotters.use_science_style()
NuEScat_dir = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_3_1_fullsample_few_recocut_recoslc'
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

def key_setx_lim(ax,key):
  """
  I realize now that we should not place cuts, but restrict x-axis (and save as seperate plot)
  """
  if key in ['true_Etheta','Etheta']:
    ax.set_xlim([-0.01,0.1])
  if key[-3:] == 'eng':
    ax.set_xlim([-0.1,4])
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
        labels.append(f'{name} ({len(xs[i]):,})')
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
        key_setx_lim(ax,key_x)
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
  for i,name in enumerate(event_names):
    x = df.loc[df['evt_type']==name,key_x].values
    #Mask certain values
    key_mask_values(x,key_x)
    x = helpers.remove_dummy_values(x,
                                    dummy_val_list=[-9999,-999,-5,9999,999],
                                    is_df=False)
    if len(x) == 0: continue #Skip event types that have no events
    xs.append(x)
    labels.append(f'{name} ({len(xs[i]):,})')
  fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)

  ax.hist(xs,label=labels,edgecolor=colors[i],stacked=True,**kwargs)
  ax.set_xlabel(key_x)
  ax.set_ylabel('Count')
  if title is not None:
    ax.set_title(title)
  ax.legend()
  plotters.set_style(ax)
  return fig,ax

def plot_masked_vals(x,masks,labels,use_kde=True,y=None,title=None,**kwargs):
  """
  Plots scatter and histograms using masked values and labels them.
  x is the data array or whateva
  y is optional and will make this a scatter plot

  """
  xs = []

  fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)
  for i,mask in enumerate(masks):
    label = f'{labels[i]} ({len(x[mask])})'
    if y is not None:
      ax.scatter(x[mask],y[mask],label=label,**kwargs)
    elif use_kde:
      ax = sns.kdeplot(x=x[mask],ax=ax,label=label,warn_singular=False,**kwargs)
  if title is not None:
    ax.set_title(title)
  plotters.set_style(ax)
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

print('-----Cleaning data')
#Mask event type and get their names
events = set_event_type(events)
events.loc[:,'nreco'] = events.nshw + events.ntrk
#events = cuts.apply_cuts(events)
event_names = categories

print('-----Seperate data types')
#Isolate data from different reconstructed objects, shws trks
lshw_keys = [key for key in events.keys() if key[:5] == 'lshw.']
lshw_keys.extend(['evt_type'])
slshw_keys = [key for key in events.keys() if key[:6] == 'slshw.']
slshw_keys.extend(['evt_type'])
ltrk_keys = [key for key in events.keys() if key[:5] == 'ltrk.']
ltrk_keys.extend(['evt_type'])
sltrk_keys = [key for key in events.keys() if key[:6] == 'sltrk.']
sltrk_keys.extend(['evt_type'])

#Get all events from dataframe
lshw = events.loc[:,lshw_keys]
slshw = events.loc[:,slshw_keys]
ltrk = events.loc[:,ltrk_keys]
sltrk = events.loc[:,sltrk_keys]

#Remove events filled with dummy values
lshw = helpers.remove_dummy_values(lshw)
slshw = helpers.remove_dummy_values(slshw)
ltrk = helpers.remove_dummy_values(ltrk)
sltrk = helpers.remove_dummy_values(sltrk)

other_keys = []
for key in events.keys():
  if key not in slshw_keys and key not in sltrk_keys and key not in ltrk_keys and key not in lshw_keys:
    if key != 'run' and key != 'subrun' and key != 'evt':
      other_keys.append(key)
other_keys.extend(['evt_type'])
other_events = events.loc[:,other_keys]
other_events = helpers.remove_dummy_values(other_events)

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
eng_bins = np.arange(0,3,0.05)
theta_bins = np.arange(0,np.pi/10,np.pi/200)
x_bins = np.arange(-250,250,10)
z_bins = np.arange(-50,550,10)
etheta_bins = np.arange(0,stop=0.02,step=0.0005)
err_bins = np.arange(-1,1,0.05)
theta_err_bins = np.arange(-2,3,0.1)
de_dx_bins = np.arange(0,10,0.1)
dens_bins = np.arange(0,50,5)
len_bins = np.arange(0,500,10)
bins_other = [int_bins,int_bins,eng_bins,theta_bins,x_bins,x_bins,z_bins,
              etheta_bins,int_bins,int_bins]

print('-----Test plots')
#plot_df_keys(other_events,keys_x=['reco_eng'])
# plot_df_keys(lshw,keys_x=['lshw.len','lshw.cnvgap','lshw.dens'])
fig,ax = plot_hist(events,'reco_eng',alpha=0.5,linewidth=2,bins=eng_bins)
plotters.save_plot(f'hist_Etheta_zoom',folder_name=folder_name,fig=fig)

print('-----Starting plots')

#Make histograms
for i,key in enumerate(events):
  if key == 'evt_type': continue

  if key[-4:] in ['dedx']:
    fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=de_dx_bins)
  elif key in ['nshw','nreco','ntrk','nstub','truenshw','truentrk','nele']:
    bins = np.arange(0,max(events.loc[:,key]),1) #arange bins into integer numbers
    fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=bins)
  elif key[-3:] == 'eng':
    fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=eng_bins)
  elif key[-4:] in ['vgap','dens']:
    fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=dens_bins)
  elif key[-3:] in ['len']:
    fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=len_bins)
  else:
    fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=20)
  if key != 'Etheta':
    key_setx_lim(ax,key)
  key = key.replace('.','_')
  plotters.save_plot(f'hist_{key}',folder_name=folder_name,fig=fig)
  if key[-6:] == 'Etheta':
    fig,ax = plot_hist(events,key,alpha=0.5,linewidth=2,bins=etheta_bins)
    key = f'{key}_zoom'
    plotters.save_plot(f'hist_{key}',folder_name=folder_name,fig=fig)
  ax.set_yscale('log')
  plotters.save_plot(f'hist_{key}_logy',folder_name=folder_name,fig=fig)
  plt.close('all')
  print(key)

plot_df_keys(other_events)
scat_keys_other = ['reco_eng','reco_theta']
plot_df_keys(other_events,keys_x=scat_keys_other,plot_kde=False,plot_scat=True)
plot_df_keys(lshw)
plot_df_keys(slshw)
plot_df_keys(ltrk)
plot_df_keys(sltrk)



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
