import os
from datetime import date
import matplotlib.pyplot as plt
import random

day = date.today().strftime("%Y_%m_%d")

def set_style(ax,legend_size=16,axis_size=16,title_size=20,tick_size=16):
  ax.tick_params(axis='x', labelsize=tick_size)
  ax.tick_params(axis='y', labelsize=tick_size)
  ax.xaxis.label.set_size(axis_size)
  ax.yaxis.label.set_size(axis_size)
  ax.title.set_size(title_size)
  ax.legend(fontsize=legend_size)

def save_plot(fname,fig=None,ftype='.png',dpi=400,folder_name=f'Plots_{day}'):
    os.system(f'mkdir -p {folder_name}')
    if fig == None:
      plt.savefig(f'{fname}{ftype}',bbox_inches = "tight",dpi=dpi)
    else:
      fig.savefig(f'{fname}{ftype}',bbox_inches = "tight",dpi=dpi)
    os.system("mv " + fname + ftype + f' {folder_name}/')

#Return dataframe from set of keys with run info as index
def get_df(tree,keys,hdrkeys=['run','subrun','evt']):
  """Input tree from uproot and keys, return df with indeces
  If dfs are different sizes, we'll return a list of dfs
  """
  copy_keys = keys.copy() #Avoid modifying input
  # if hdrkeys not in copy_keys:
  #   copy_keys.extend(hdrkeys)
  df = tree.arrays(copy_keys,library='pd')
  if isinstance(df,tuple): #If it's a tuple, we'll rename each df, and return the list of them
    dfs = []
    for tempdf in df:
      tempdf.set_index(hdrkeys)
      tempdf = tempdf.sort_index()
      dfs.append(tempdf)
    return dfs
  else: #Returns single df
    df.set_index(hdrkeys)
    df = df.sort_index()
    return df

def get_random_colors(n=10):
  # Generate n random colors
  return [tuple(random.random() for _ in range(3)) for _ in range(n)]