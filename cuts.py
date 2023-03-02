"""
Preceed each cut with a c
"""


def cEtheta2(e_theta,cut_val=0.004):
  #GeV^2 rad
  return e_theta<cut_val

def cERazzle(erazzle,cut_val=0.6):
  return erazzle<cut_val

def cNShw(nshw,cut_val=1):
  return nshw == cut_val

def cNTrk(ntrk,cut_val=0):
  return ntrk == cut_val
    
def cTheta(theta,cut_val=60e-3):
  #rad
  return theta<cut_val

def cNEle(nele,cut_val=1):
  return nele==cut_val

def apply_cuts(events):
  """
  Apply cuts to dataframe
  """
  for ind,row in events.iterrows():
    #Cuts return true if they're passed, tweak these to modify values
    events.loc[ind,'cEtheta2_true'] = cEtheta2(row.true_Etheta)
    events.loc[ind,'cNEle_true'] = cNEle(row.nele)
    events.loc[ind,'cNShw_true'] = cNShw(row.truenshw)
    events.loc[ind,'cNTrk_true'] = cNTrk(row.truentrk)
    events.loc[ind,'cNShw_reco'] = cNShw(row.nshw)
    events.loc[ind,'cNTrk_reco'] = cNTrk(row.ntrk)
    events.loc[ind,'cERazzle_lshw'] = cERazzle(row['lshw.electron'])
    events.loc[ind,'cEtheta2_reco'] = cEtheta2(row['Etheta'])
  return events