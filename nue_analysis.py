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
NuEScat_dir = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_2_27_fullsample_nue_no_cuts_truthslc'
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

def plot_all_keys(df,ftype='.png'):
  """
  plot_all_keys(df)

  Parameters:
  df (pandas dataframe): dataframe that contains the data to be plotted

  Returns:
  None
  """
  counter = 0
  fnames = []
  alphas = [1,0.4,0.4,0.4,0.4,0.4]
  os.system(f'mkdir -p {folder_name}')
  for key_x in df.keys():
    for key_y in df.keys():
      #if key_x != key_y or key_x == 'evt_type': continue
      if key_y == 'evt_type': continue
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
            alphas[i] = 0.05
        x,drop_ind = helpers.remove_dummy_values(df.loc[df['evt_type']==name,key_x].values,
                                        dummy_val_list=[-9999,-999,-5,9999,999],
                                        is_df=False,
                                        return_drop_indeces=True)
        y = np.delete(df.loc[df['evt_type']==name,key_y].values,drop_ind) #use indeces from x for same size scatter plot
        xs.append(x)
        ys.append(y)
        labels.append(f'{name} ({len(xs[i])})')
        if key_x != key_y:
          sns.scatterplot(x=xs[i],y=ys[i],label=labels[i],ax=ax,alpha=alphas[i])
          ax.set_xlabel(f'{key_x}')
          ax.set_ylabel(f'{key_y}')
          if key_x == 'Etheta':
            ax.set_xlim([0,0.02])
          if key_x[-4:] == 'dedx':
            ax.set_xlim([0,10])
          prefix = 'scat'
        if key_x == key_y and key_x != 'evt_type':
          #weights = [abs(x/len(x)) for x in xs] #Weight each bin to its weight. Each histogram bin should be one
          #ax.hist(weights,label=labels,alpha=0.3,edgecolor=colors[i],stacked=True,density=False)#,weights=weights)
          #sns.histplot(x=xs,ax=ax,label=labels,alpha=alphas,common_norm=False,
          #common_bins=False,bins=10,color=colors)
          ax = sns.kdeplot(x=xs[i],ax=ax,label=labels[i],warn_singular=False,alpha=alphas[i])
          ax.set_xlabel(f'{key_x}')
          ax.set_ylabel('Counts')
          prefix = 'dist'
        
        
      ax.legend()
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
    x = helpers.remove_dummy_values(df.loc[df['evt_type']==name,key_x].values,
                                    dummy_val_list=[-9999,-999,-5,9999,999],
                                    is_df=False)
    #if len(x) == 0: continue #Skip event types that have no events
    xs.append(x)
    labels.append(f'{name} ({len(xs[i])})')
  fig,ax = plt.subplots(figsize=(6,6),tight_layout=True)

  ax.hist(xs,label=labels,edgecolor=colors[i],stacked=True,**kwargs)
  ax.set_xlabel(key_x)
  ax.set_ylabel('Count')
  if title is not None:
    ax.set_title(title)
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
  tree = uproot.open(f'{NuEScat_dir}/{fname}:rectree;1')
  #use numpy library to convert numpy array to dataframe
  df = helpers.get_df(tree,
                      tree.keys(),
                      hdrkeys=['run','subrun','evt'],
                      library='np')
  dfs.append(df)
events = pd.concat(dfs,axis=1)

#Mask event type and get their names
events = set_event_type(events)
events.loc[:,'nreco'] = events.nshw + events.ntrk
events = cuts.apply_cuts(events)
#events.to_csv(f'{NuEScat_dir}/events_raw.csv')
#events = pd.read_csv(f'{NuEScat_dir}/events_raw.csv')
#events = events[(events != -9999).all(1)] #Remove dummy values
#events = events[(events != -999).all(1)] #Remove dummy values
event_names = categories

#Isolate data from different reconstructed objects, shws trks
lshw_keys = [key for key in events.keys() if key[:5] == 'lshw.']
lshw_keys.extend(['evt_type'])
slshw_keys = [key for key in events.keys() if key[:6] == 'slshw.']
slshw_keys.extend(['evt_type'])

lshw = events.loc[:,lshw_keys]
slshw = events.loc[:,slshw_keys]

lshw = helpers.remove_dummy_values(lshw)
#lshw = lshw.loc[lshw.loc[:,'lshw.dedx']<100] #Remove all huge dedx values

slshw = helpers.remove_dummy_values(slshw)
#slshw = slshw.loc[slshw.loc[:,'slshw.dedx']<100] #Remove all huge dedx values

ltrk_keys = [key for key in events.keys() if key[:5] == 'ltrk.']
ltrk_keys.extend(['evt_type'])
sltrk_keys = [key for key in events.keys() if key[:6] == 'sltrk.']
sltrk_keys.extend(['evt_type'])

ltrk = events.loc[:,ltrk_keys]
sltrk = events.loc[:,sltrk_keys]

ltrk = helpers.remove_dummy_values(ltrk)
#ltrk = ltrk.loc[ltrk.loc[:,'ltrk.dedx']<100] #Remove all huge dedx values

sltrk = helpers.remove_dummy_values(sltrk)
#sltrk = sltrk.loc[sltrk.loc[:,'sltrk.dedx']<100] #Remove all huge dedx values

other_keys = []
for key in events.keys():
  if key not in slshw_keys and key not in sltrk_keys and key not in ltrk_keys and key not in lshw_keys:
    if key != 'run' and key != 'subrun' and key != 'evt':
      other_keys.append(key)
other_keys.extend(['evt_type'])
other_events = events.loc[:,other_keys]
other_events = helpers.remove_dummy_values(other_events)

plot_all_keys(lshw)
plot_all_keys(other_events)

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
eng_bins = np.arange(0,1.5,0.05)
theta_bins = np.arange(0,np.pi/10,np.pi/200)
x_bins = np.arange(-250,250,10)
z_bins = np.arange(-50,550,10)
etheta_bins = np.arange(0,stop=0.01,step=0.0001)
err_bins = np.arange(-1,1,0.05)
theta_err_bins = np.arange(-2,3,0.1)
de_dx_bins = np.arange(0,5,0.1)
bins_other = [int_bins,int_bins,eng_bins,theta_bins,x_bins,x_bins,z_bins,
              etheta_bins,int_bins,int_bins]

#Make plots
for rl,tl,ra,ta,bins in zip(reco_labels_other,true_labels_other,reco_array_other,true_array_other,bins_other):
  if rl in ['nshw','ntrk','nreco']:
    fig,ax = hist2d_true_reco(ra,ta,rl,tl,title=r'$\nu+e$',bins=bins,cmap='viridis_r',label_boxes=True)
    plotters.save_plot(f'hist2d_{tl}_{rl}',folder_name=folder_name,fig=fig)
  # elif rl == 'reco_eng':  
  #   fig,ax = hist2d_true_reco(ra,ta,rl,tl,title=r'$\nu+e$',plot_line=True,bins=bins)
  #   #ax.plot([0,1],[0,2],ls='--',color='g')
  #   #ax.plot([0,2/3],[0,2],ls='--',color='c')
  #   plotters.save_plot(f'hist2d_{tl}_{rl}',folder_name=folder_name,fig=fig)
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

for i,key in enumerate(lshw_keys):
  if key == 'lshw.dedx' or key == 'lshw.cnvgap':
    fig,ax = plot_hist(lshw,key,alpha=0.5,linewidth=2,bins=de_dx_bins)
  else:
    fig,ax = plot_hist(lshw,key,alpha=0.5,linewidth=2,bins=20)
  key = key.replace('.','_')
  plotters.save_plot(f'hist_{key}',folder_name=folder_name,fig=fig)
  plt.close('all')

for i,key in enumerate(other_keys):
  fig,ax = plot_hist(other_events,key,alpha=0.5,linewidth=2,bins=20)
  key = key.replace('.','_')
  plotters.save_plot(f'hist_{key}',folder_name=folder_name,fig=fig)
  ax.set_yscale('log')
  plotters.save_plot(f'hist_{key}_logy',folder_name=folder_name,fig=fig)
  plt.close('all')




#Plot masked values
# lshw_inds = lshw.index
# passed = events.cNShw_reco[lshw_inds] & events.cNTrk_reco[lshw_inds] & events.cERazzle_lshw[lshw_inds] & events.cEtheta2_reco[lshw_inds] 
# mask_cuts = {'Passed':passed,
#              'Failed':[not ele for ele in passed]}
# masks = list(mask_cuts.values())
# labels = list(mask_cuts.keys())
# cnt = 0
# for i,xkey in enumerate(lshw_keys):
#   for j,ykey in enumerate(lshw_keys):
#     if xkey == 'evt_type' or ykey == 'evt_type': continue
#     x = lshw.loc[:,xkey].values
#     y = lshw.loc[:,ykey].values
#     if xkey == ykey:
#       fig,ax = plot_masked_vals(x,masks,labels,title=r'$\nu+e$')
#       ax.set_xlabel(f'{xkey}')
#       plotters.save_plot(f"maskhist_{xkey.replace('.','_')}",folder_name=folder_name,fig=fig)
#     else: 
#       fig,ax = plot_masked_vals(x,masks,labels,y=y,title=r'$\nu+e$',alpha=0.3)
#       ax.set_xlabel(f'{xkey}')
#       ax.set_ylabel(f'{ykey}')
#       plotters.save_plot(f"maskscat_{xkey.replace('.','_')}_{ykey.replace('.','_')}",folder_name=folder_name,fig=fig)
#     plt.close('all')
    
#print(other_keys)
#plot_all_keys(other_events)
#plot_all_keys(lshw)
plot_all_keys(slshw)
plot_all_keys(ltrk)
plot_all_keys(sltrk)

