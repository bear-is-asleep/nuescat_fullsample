from datetime import date
import random
import numpy as np
import pandas as pd

day = date.today().strftime("%Y_%m_%d")

#Return dataframe from set of keys with run info as index
def get_df(tree,keys,hdrkeys=['run','subrun','evt'],library='pd'):
  """Input tree from uproot and keys, return df with indeces
  If dfs are different sizes, we'll return a list of dfs
  """
  copy_keys = keys.copy() #Avoid modifying input
  # if hdrkeys not in copy_keys:
  #   copy_keys.extend(hdrkeys)
  df = tree.arrays(copy_keys,library=library)
  if library == 'pd':
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
  elif library == 'np':
    ret_df = pd.DataFrame() #Make empty dataframe
    for key,item in df.items(): #numpy library returns a dict of arrays, make df using this
      ret_df[key] = item[0]
    ret_df = ret_df.set_index(hdrkeys)
    ret_df = ret_df.sort_index()
    return ret_df #still return a dataframe (im persistent)


def get_random_colors(n=10):
  # Generate n random colors
  return [tuple(random.random() for _ in range(3)) for _ in range(n)]

def remove_dummy_values(events,dummy_val_list=[-9999,-999,999,9999],is_df=True,
                        return_drop_indeces=False):
  """
  Removes values in df or array from dummy list
  is_df = True for dataframe, false for numpy array (must be 1D)
  """
  if is_df:
    #Remove dummy values from events df
    for val in dummy_val_list:
      events = events[(events != val).all(1)]
  else: 
    drop_ind = []
    for val in dummy_val_list:
      drop_condition = np.where(events==val) #Drop all values that are the dummy value
      if len(drop_condition[0]) == 0: continue #skip if we don't find any values to drop
      drop_ind.extend(list(drop_condition[0]))
      events = np.delete(events,drop_condition)
  if return_drop_indeces:
    return events,drop_ind
  else:
    return events

def correct_angles(events):
  """
  Return corrected angles based on vertex location
  """
def confine_points(x):
  return 1
