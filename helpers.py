import os
from datetime import date
import matplotlib.pyplot as plt
import random

day = date.today().strftime("%Y_%m_%d")

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
      tempdf = tempdf.set_index(hdrkeys)
      tempdf = tempdf.sort_index()
      dfs.append(tempdf)
    return dfs
  else: #Returns single df
    df = df.set_index(hdrkeys)
    df = df.sort_index()
    return df

def get_random_colors(n=10):
  # Generate n random colors
  return [tuple(random.random() for _ in range(3)) for _ in range(n)]

def remove_dummy_values(events,dummy_val_list=[-9999,-999,-5]):
  #Remove dummy values from events df
  for val in dummy_val_list:
    events = events[(events != val).all(1)]
  return events
