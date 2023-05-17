from datetime import date
import random
import numpy as np
import pandas as pd

day = date.today().strftime("%Y_%m_%d")

colors = ['blue','purple','yellow','red','green','brown','orangered','black']
categories = [r'$\nu + e$',r'NC $\pi^0$','NC',r'CC $\nu_\mu$',r'CC $\nu_e$','Dirt','Cosmics','Other']

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
      events = events[(events != val)]
  else: 
    drop_ind = []
    for val in dummy_val_list:
      drop_condition = np.where(events==val) #Drop all values that are the dummy value
      if len(drop_condition[0]) == 0: continue #skip if we don't find any values to drop
      drop_ind.extend(list(drop_condition[0]))
    events = np.delete(events,drop_ind)
  if return_drop_indeces:
    return events,drop_ind
  else:
    return events

def set_event_type(events):
  """
  set_event_type(events)

  This function takes a dataframe (events) as input, and sets the 
  'evt_type' column to a string value based on the numerical value of
  'evt_type' column. The numerical values are mapped to string values: 
  0: 'NuEScat', 1: 'NCPi0', 2: 'NC', 3: 'CCNuMu', 4: 'CCNuE', 5: 'Dirt', 6: 'Cosmic', 7: 'Other'.

  Parameters:
  events (pandas dataframe): dataframe that contains the event data

  Returns:
  events (pandas dataframe): modified dataframe with the 'evt_type' column as string values
  """
  
  for i,cat in enumerate(categories):
    events.loc[events.loc[:,'evt_type'] == i,'evt_type'] = cat
  return events

def calc_etheta_reco(events,prefix):
  """
  Calc etheta for reco object, and mask out bad values
  """
  events.loc[:,f'{prefix}Etheta'] = events.loc[:,f'{prefix}eng']* events.loc[:,f'{prefix}angle']**2
  mask = abs(events.loc[:,f'{prefix}Etheta']) > 9999 #Mask these values to -9999
  mask_inds = events.index[mask]
  events.loc[mask_inds,f'{prefix}Etheta'] = -9999 #Mask to -9999
  return events

def combine_bins_sum(counts,edges,iterations=1):
  """
  Takes counts and edges and combines adjacent bins
  iterations is how many times to combine the bins
  """
  for _ in range(iterations):
    new_edges = np.zeros(int(len(edges) / 2 + 1)) # Create an empty array to hold new bin edges
    new_counts = np.zeros(int(len(counts)/2)+len(counts)%2) # Create an empty array to hold new bin counts
    for i in range(len(edges)):
        ind = int(i/2)
        if i%2 == len(edges)%2:
          new_counts[ind] = counts[i-1] + counts[i]
        else:
          new_edges[ind] = edges[i]
    new_edges[-1] = edges[-1] # Add the last bin count
    edges = new_edges
    counts = new_counts
  return counts,edges

def frac_syst_stat_bias_heights(values,edges,weights):
  """
  reweight flux for each energy bin by weight

  universes should be binned by energy

  Returns constrained and unconstrained bin heights, and fractional uncertainties
  """
  #Make arrays to fill
  height_constrained = np.zeros(len(edges)-1)
  height = height_constrained.copy()

  #Resolution
  fracunc_constrained = height_constrained.copy()
  fracunc = height_constrained.copy()

  #Bias
  fracbias_constrained = height_constrained.copy()
  fracbias = height_constrained.copy()

  #Stat
  fracstat_constrained = height_constrained.copy()
  fracstat = height_constrained.copy()

  for i in range(height_constrained.shape[0]):
    vals = values[1:,i] #Grab universe bin heights within energy bin
    nom = values[0,i] #Central value flux
    #Get heights
    height_constrained[i] = np.average(vals, weights=weights)
    height[i] = np.average(vals)
    #Get resolutions
    fracunc_constrained[i] = np.sqrt(np.average((vals - height_constrained[i]) ** 2, weights=weights))/height_constrained[i]
    fracunc[i] = np.sqrt(np.average((vals - height[i]) ** 2))/height[i]
    #Get biases
    fracbias_constrained[i] = np.sqrt(np.average((vals - nom) ** 2,weights=weights))/nom
    fracbias[i] = np.sqrt(np.average((vals - nom) ** 2))/nom
    #Get stat
    fracstat_constrained[i] = np.average(np.sqrt(vals),weights=weights)/height[i]
    fracstat[i] = np.average(np.sqrt(vals))/height[i]
  return height_constrained,height,fracunc_constrained,fracunc,fracbias_constrained,fracbias,fracstat_constrained,fracstat

def drop_matching_indices(lst, indices_to_check):
  """
  Takes a list and a list of indices to check, and returns a new list with the values
  that do not match those indices.
  """
  return [value for index, value in enumerate(lst) if index not in indices_to_check]

def Q2_nuecc(E_e,theta_e):
  """
  Q2 = 2mn(Ev-Ee)
  Ev = (mnEe-me^2/2)/(mn-Ee+Pe cos0_e)
  """
  m_e = 0.511e-3 #GeV
  m_n = 0.938 #GeV
  p_e = np.sqrt(E_e**2-m_e**2)
  E_v = (m_n*E_e-m_e**2/2)/(m_n-E_e+p_e*np.cos(theta_e))
  return 2*m_n*(E_v-E_e)