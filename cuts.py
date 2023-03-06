"""
Preceed each cut with a c
"""

import numpy as np

def efficiency_purity_f1_calc(events,mask):
  """
  Calc efficiency for cuts, apply a mask to get lost values
  """
  int_events = len(events)
  intsig_events = np.count_nonzero(events.evt_type == r'$\nu + e$')
  events = events[mask] #mask to get only events that pass cut
  sig_events = np.count_nonzero(events.evt_type == r'$\nu + e$')
  tot_events = len(events)

  #Calc eff,pur,f1
  eff = sig_events/intsig_events
  if tot_events == 0:
    pur = 0
    f1 = 1/2*eff
  else:
    pur = sig_events/tot_events
    f1 = 1/2*(1/pur+1/eff) #f1 score
  return pur,eff,f1

def cEtheta2(events,cut_key,bool_key='cEtheta2',cut_val=0.01):
  """
  Return bool_key which indicates whether the event passed the cut
  """
  #Use abs since some mask values are really negative :/
  events[bool_key] = abs(events.loc[:,cut_key])<cut_val
  return events

def cLShwLen(events,cut_key,bool_key='cRecoLen',cut_val=1.):
  """
  Return bool_key which indicates whether the event passed the cut
  """
  #Use abs since some mask values are really negative :/
  events[bool_key] = abs(events.loc[:,cut_key])<cut_val
  return events

def select_electron(reco_objs,priors=[1,0,0,0]):
  """
  Select electron using object with lowest theta

  reco_objs is a list of dataframes containing info about reconstructed objects
  priors contain prior probability of being the electron
  """
  #In a study with all of the reco nu+e, the shower is the electron 143/157 times
  return 0 #assume order is lshw,lrtk,slshw,slrtk
def electron_etheta_theta_eng(reco_objs,distance_threshold=5,angle_theshold=0.1):
  """
  Take in all reco objects, find electron, and return the angle, energy and etheta
  distance threshold so objects are within some reasonable distance between the two
  angle threshold so they travel in appx. the same direction
  """
  eid = select_electron(reco_objs) #electron id
  electron = reco_objs[eid] #electron object

  #This assumes everything is selected as a track, you'd probably need to update if not assuming this
  electron_start_keys = [key for key in electron.keys() if 'start' in key] #get start keys
  electron_angle_key = 'lshw.angle'
  electron_start = electron.loc[:,electron_start_keys]
  for i,reco_obj in enumerate(reco_objs):
    if i == eid: continue#dont double count
    #get keys
    start_keys = [key for key in reco_obj.keys() if 'start' in key] #get start keys
    angle_key = [key for key in reco_obj.keys() if key[-5:] == 'angle' and 'true' not in key] #get reco angle
    eng_key = [key for key in reco_obj.keys() if key[-3:] == 'eng' and 'true' not in key] #get reco angle
    
    #get values
    reco_start = reco_obj.loc[:,start_keys]
    distance = np.linalg.norm(electron_start-reco_start,axis=1)
    #angle_difference = 



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

