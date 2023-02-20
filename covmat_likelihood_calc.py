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
state_dir = f'/sbnd/data/users/brindenc/analyze_sbnd/nue/states/2022A/2023_1_30_systs_5bins'

fname = 'state_all.root' #To get all universes and nominal one
hist = uproot.open(f'{state_dir}/{fname}') #Open histograms

#GEt bin values and edges
keys = hist.keys()
edges = hist[keys[0]].axis().edges() #They all have the same bin widths
values = np.zeros((len(keys),len(edges)-1)) #One less bin edge shape
for i,key in enumerate(keys): #Iterate over all universes
  values[i] = hist[key].values()

#Extract nominal and universe values
nom = values[0] #I temporarily added the :-1 because the last bin is always zero
universes = values[1:]
n_universes = universes.shape[0]

print(universes)
plt.hist(np.sum(universes,axis=1))
plt.savefig('tests/universe.png')

#Calcu cov mat using Vij = 1/N sum_i((xu_i-xn_i)(xu_j-xn_j))
covmat = np.zeros((len(nom),len(nom)))

for _,uni in enumerate(universes):
  for i,_ in enumerate(nom):
    for j,_ in enumerate(nom):
      covmat[i][j] += (uni[i]-nom[i])*(uni[j]-nom[j])

covmat/=n_universes

# Display the covmat 
if True:
  fig,ax = plt.subplots(figsize=(8,8),tight_layout=True)
  plt.imshow(covmat,cmap='seismic')
  # add text to display the value of each bin
  for i in range(covmat.shape[0]):
    for j in range(covmat.shape[1]):
      plt.text(j, i, f'{covmat[i,j]:.1f}', ha='center', va='center', color='white')
  plt.xticks(np.arange(len(edges)-1),edges[:-1]) #Label ticks with energy spectrum
  plt.yticks(np.arange(len(edges)-1),edges[:-1]) #Label ticks with energy spectrum
  plt.colorbar()
  helpers.save_plot('covmat',folder_name=state_dir)
print('covmat made')

#Calc likelihood W = P(N_{\nu+e}|M) = \frac{1}{(2\pi)^{K/2}}\frac{1}{\sqrt{\Sigma_N}}
# \textrm{exp}\left( -\frac12 (\textbf{N}-\textbf{M})^T\Sigma_N^{-1} (\textbf{N}-\textbf{M})\right)
#χ2 =(N−M)T Σ−1(N−M)
print(covmat)
dof = len(values[0]) #degrees of freedom
chi_squared_arr = np.zeros(n_universes)
normalization = 1/((2*np.pi)**(dof/2)*np.linalg.det(covmat)**(1/2))
W_arr = np.zeros(n_universes)
for i,uni in enumerate(universes):
  chi_squared_arr[i] = (nom-uni).T@np.linalg.inv(covmat)@(nom-uni)
  W_arr[i] = normalization*np.exp(-0.5*chi_squared_arr[i])

#print(W_arr,chi_squared_arr)
if True:
  plt.hist(W_arr)
  plt.savefig('tests/ws.png')

  plt.hist(chi_squared_arr)
  plt.savefig('tests/x2.png')
