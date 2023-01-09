import sys
sys.path.append('/sbnd/app/users/brindenc/mypython') #My utils path
from bc_utils.CAFana import pic as CAFpic
from bc_utils.CAFana import plotters as CAFplotters
from bc_utils.utils import pic,plotters
from time import time
import numpy as np
import pandas as pd
import uproot
import ROOT

plt.style.use(['science','no-latex'])
day = date.today().strftime("i%Y_%m_%d")

#Params import
from CAFdata import * #Its okay to do this, they're all just variable names anyway

