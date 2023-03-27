"""
Preceed each cut with a c
"""

import numpy as np
from time import time
import logging
import helpers

def efficiency_purity_f1_calc(events,mask,intsig_events=None,tot_events=None,
                              max_pur=1,max_eff=1):
  """
  Calc efficiency for cuts, apply a mask to get lost values
  """
  if intsig_events is None:
    intsig_events = np.count_nonzero(events.evt_type == r'$\nu + e$')
  events = events[mask] #mask to get only events that pass cut
  sig_events = np.count_nonzero(events.evt_type == r'$\nu + e$')
  if tot_events is None:
    tot_events = len(events)

  #Calc eff,pur,f1
  
  eff = sig_events/intsig_events
  if tot_events == 0 and sig_events == 0:
    return 0,0,0,sig_events,tot_events
  elif sig_events == 0:
    return 0,0,0,sig_events,tot_events
  elif tot_events == 0:
    pur = 0
    f1 = 1/(2*eff)
  else:
    pur = sig_events/tot_events
    f1 = 1/2*(1/pur+1/eff) #f1 score
  return pur,eff,f1,sig_events,tot_events

def f1_calc(pur,eff,max_pur=1,max_eff=1):
  """
  Calcs f1 weighted by max purity and efficiency
  """
  return 1/2*(max_pur/pur+max_eff/eff)

def cBoxCut(events,cut_key,cut_val,bool_key,equality='lessthan',
  use_abs=True):
  """
  Return bool_key which indicates whether the event passed the cut
  equality determines if we use less than or greater than cut val
  cut_val is a list of cut_vals for multiple values
  --lessthan : keep events less than
  --greaterthan : keep events greater than
  --isbetween : keep events less than cut_val[0], greater than cut_val[1]
  use_abs : use absolute value to cut on
  """
  events.loc[:,bool_key] = True #mask as true for now
  if use_abs:
    vals = abs(events.loc[:,cut_key])
  else:
    vals = events.loc[:,cut_key]
  vals = helpers.remove_dummy_values(vals,
                                       dummy_val_list=[-9999,-999,999,9999]) #Remove dummy vals for search, don't cut on masked out values
  mask_inds = vals.index
  if equality == 'lessthan':
    events.loc[mask_inds,bool_key] = vals<cut_val
  elif equality == 'greaterthan':
    events.loc[mask_inds,bool_key] = vals>cut_val
  elif equality == 'isbetween':
    events.loc[mask_inds,bool_key] = vals<cut_val[1] and vals>cut_val[0]
  else:
    logging.warning(f'{equality} is not a defined equality')
  return events

def cNRecoCut(events,cut_key,nreco_key,nreco,cut_val,bool_key,equality='lessthan',nreco_equality='equal',
              use_abs=True):
  """
  Returns events with bool key. Makes cuts on different values of ntrk,nstub,nshw, etc
  by accessing cut_dict
  nreco number of reco objects to consider
  cut_val value to cut on
  last_key is last key to check, at this one we just check for values greater than nreco
  nreco_equality asks if we want to include all nreco equal to, lessthan or greater than the 
  value nreoc specifies
  """
  events.loc[:,bool_key] = True #mask as true for now
  if nreco_equality == 'equal':
    mask = events.loc[:,nreco_key] == nreco #Mask to events with this value of nreco
  elif nreco_equality == 'lessthan':
    mask = events.loc[:,nreco_key] < nreco #Mask to events with this value of nreco
  elif nreco_equality == 'greaterthan':
    mask = events.loc[:,nreco_key] > nreco #Mask to events with this value of nreco
  if use_abs:
    vals = abs(events[mask].loc[:,cut_key])
  else:
    vals = events[mask].loc[:,cut_key]
  vals = helpers.remove_dummy_values(vals,
                                       dummy_val_list=[-9999,-999,5,999,9999]) #Remove dummy vals for search, don't cut on masked out values
  mask_inds = vals.index
  #logging.info(f'vals = {vals}')
  if equality == 'lessthan':
    events.loc[mask_inds,bool_key] = vals<cut_val
  elif equality == 'greaterthan':
    events.loc[mask_inds,bool_key] = vals>cut_val
  elif equality == 'isbetween':
    events.loc[mask_inds,bool_key] = vals<cut_val[1] and vals>cut_val[0]
  else:
    logging.warning(f'{equality} is not a defined equality')
  #events.loc[:,bool_key] = events.loc[:,bool_key].fillna(True).values #Give values not equal to the nreco we're looking for a pass
  #logging.info(f'{events.loc[:,bool_key].head(30)}')
  return events
    


def mask_cosmics(events,time_threshold=None,angle_threshold=None,prefixes=['ltrk.','sltrk.'],
                 dummy_val=-9999.):
  """
  Muons that have crt interactions are muons from outside the detector. 
  We should mask these and make cuts on the them as if the cosmics dont exist

  time_threshold/angle_threshold should be a vector of greater than less than values
  
  """

  for _,prefix in enumerate(prefixes):
    prefix_keys = [key for key in events.keys() if prefix == key[:len(prefix)]] #Get all keys with prefix in front
    
    angle_mask = events.loc[:,f'{prefix}crttrk.angle'].values != dummy_val #if there is a crttrk this will be true
    time_mask = events.loc[:,f'{prefix}crttrk.time'].values != dummy_val #if there is a crttrk this will be true
    if time_threshold is not None: #constrain to values less than max, greater than min
      upper_bound = events.loc[:,f'{prefix}crttrk.time'].values > time_threshold[1]
      lower_bound = events.loc[:,f'{prefix}crttrk.time'].values < time_threshold[0]
      
      time_mask = np.logical_and(time_mask,np.logical_or(upper_bound,lower_bound)) #Combine lower and upper bound constraints
    if angle_threshold is not None: #constrain to values less than max, greater than min
      upper_bound = events.loc[:,f'{prefix}crttrk.angle'].values > angle_threshold[1]
      lower_bound = events.loc[:,f'{prefix}crttrk.angle'].values < angle_threshold[0]

      angle_mask = np.logical_and(angle_mask,np.logical_or(upper_bound,lower_bound)) #Combine lower and upper bound constraints
    mask = np.logical_and(angle_mask,time_mask) #Find all cosmics that satisfy criteria
    mask_inds = events.index.values[mask]
    events.loc[mask_inds,prefix_keys] = -9999. #set to dumby value
    events.loc[mask_inds,'ntrk'] -= 1 #subtract one from tracks
    events.loc[mask_inds,'nreco'] -= 1 #subtract one from total reco objects
  #logging.info(f'Mask cosmics ltrk -- {events.loc[:,"ltrk.crttrk.time"].head(40)}')
  #logging.info(f'Mask cosmics sltrk -- {events.loc[:,"sltrk.crttrk.time"].head(40)}')
  return events #This will set all cosmic track values to -9999, indicating it shouldn't be included in cut
def select_electron(events,prefixes=['ltrk.','slshw.','sltrk.'],distance_threshold=None,angle_threshold=None,erazzle_threshold=0.3,
                    dummy_val = -9999):
  """
  Select electron using object with lowest theta

  reco_objs is a list of dataframes containing info about reconstructed objects
  priors contain prior probability of being the electron, for now the leading shower is assumed to be the electron
  distance_threshold should be value which compares the start location of other reco objects
  angle_threshold should be value which compares the angle of the object in question to the electron
  """
  #In a study with all of the reco nu+e, the shower is the electron 143/157 times at truth, or 126/132 at reco (CRUMBS slice)
  lshw_keys = [key for key in events.keys() if key[:len(prefixes[0])] == 'lshw.'] #Get shower keys
  electron_keys = ['ele.'+key[5:] for key in lshw_keys] #Make electron keys
  events.loc[:,electron_keys] = events.loc[:,lshw_keys].values #Assume leading shower is electron

  #Mask out low erazzle showers
  #events.loc[:,'ele.electron'] < 

  #Electron vertex
  start_keys = ['ele.start.x','ele.start.y','ele.start.z']
  ele_start = events.loc[:,start_keys].values

  #Electron angle
  angle_key = 'ele.angle'
  ele_angle = events.loc[:,angle_key].values

  #Electron eng
  eng_key = 'ele.eng'
  ele_eng = events.loc[:,eng_key].values

  #Electron etheta
  lshw_ele_mask = events.loc[:,'lshw.angle'] != -9999
  lshw_ele_inds = events.index[lshw_ele_mask]

  #logging.info(f'{events.loc[:,electron_keys].head()}')

  for _,prefix in enumerate(prefixes):
    if distance_threshold is None and angle_threshold is None: break
    prefix_keys = [key for key in events.keys() if prefix == key[:len(prefix)]] #Get all keys with prefix in front

    #Get start keys for object
    reco_obj_start_keys = [prefix+key[4:] for key in start_keys]
    reco_obj_start = events.loc[:,reco_obj_start_keys].values

    #Get angle for object
    reco_obj_angle_key = prefix+'angle'
    reco_obj_angle = events.loc[:,reco_obj_angle_key].values

    #Get energy for object
    reco_obj_eng_key = prefix+'eng'
    reco_obj_eng = events.loc[:,reco_obj_eng_key].values

    #We have no choice but to loop through all of these values, to not calculate dummy_vals
    mask_inds = []
    start_time = time()
    for rs,es,ra,ea,re,ee,ind in zip(reco_obj_start,ele_start,reco_obj_angle,ele_angle,reco_obj_eng,ele_eng,events.index):
      if dummy_val in rs or dummy_val in es: continue #skip comparisons between events that don't exist
      distance = np.linalg.norm(rs-es)
      angle_difference = abs(ra-ea)
      #Conditions for being part of electron, add additional condition that the reco object has less energy than the electron
      if distance < distance_threshold and angle_difference < angle_threshold and re<ee: 
        events.loc[ind,eng_key] = re+ee #Add energy from other object
        mask_inds.append(ind)
    end_time = time()
    logging.info(f'Masking time for {prefix} in electron selection: {end_time-start_time:.1f} (s)')
  
  events.loc[:,'ele.Etheta'] = -9999 #Set dummy value
  events.loc[lshw_ele_inds,'ele.Etheta'] = events.loc[:,'ele.eng']*events.loc[:,'ele.angle']**2 #add etheta^2 info
  #logging.info(f'{events.loc[:,electron_keys].head()}')
  return events
  



# def cEtheta2(e_theta,cut_val=0.004):
#   #GeV^2 rad
#   return e_theta<cut_val

# def cERazzle(erazzle,cut_val=0.6):
#   return erazzle<cut_val

# def cNShw(nshw,cut_val=1):
#   return nshw == cut_val

# def cNTrk(ntrk,cut_val=0):
#   return ntrk == cut_val
    
# def cTheta(theta,cut_val=60e-3):
#   #rad
#   return theta<cut_val

# def cNEle(nele,cut_val=1):
#   return nele==cut_val

# def apply_cuts(events):
#   """
#   Apply cuts to dataframe
#   """
#   for ind,row in events.iterrows():
#     #Cuts return true if they're passed, tweak these to modify values
#     events.loc[ind,'cEtheta2_true'] = cEtheta2(row.true_Etheta)
#     events.loc[ind,'cNEle_true'] = cNEle(row.nele)
#     events.loc[ind,'cNShw_true'] = cNShw(row.truenshw)
#     events.loc[ind,'cNTrk_true'] = cNTrk(row.truentrk)
#     events.loc[ind,'cNShw_reco'] = cNShw(row.nshw)
#     events.loc[ind,'cNTrk_reco'] = cNTrk(row.ntrk)
#     events.loc[ind,'cERazzle_lshw'] = cERazzle(row['lshw.electron'])
#     events.loc[ind,'cEtheta2_reco'] = cEtheta2(row['Etheta'])
#   return events

