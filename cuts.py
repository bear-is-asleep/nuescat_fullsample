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

