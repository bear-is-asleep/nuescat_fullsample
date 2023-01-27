import pandas as pd
import seaborn as sns
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import sys
sys.path.append('../')
from CAFdata import *
import helpers
from datetime import date

day = date.today().strftime("%Y_%m_%d")
state_dir = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_26_systs'

fname = 'state_all.root' #To get all universes and nominal one
hist = uproot.open(f'{state_dir}/{fname}') #Open histograms

#GEt bin values and edges
keys = hist.keys()
edges = hist[keys[0]].axis().edges() #They all have the same bin widths
values = np.zeros((len(keys),len(edges)-1)) #One less bin edge shape
for i,key in enumerate(keys): #Iterate over all universes
  values[i] = hist[key].values()

#Extract nominal and universe values
nom = values[0]
universes = values[1:]
n_universes = universes.shape[0]

#Calcu cov mat using Vij = 1/N sum_i((xu_i-xn_i)(xu_j-xn_j))
covmat = np.zeros((len(nom),len(nom)))

for _,uni in enumerate(universes):
  for i,_ in enumerate(nom):
    for j,_ in enumerate(nom):
      covmat[i][j] += (uni[i]-nom[i])*(uni[j]-nom[j])

covmat/=n_universes

# Display the covmat 
plt.imshow(covmat)
plt.colorbar()
plt.savefig('tests/covmat.png')


