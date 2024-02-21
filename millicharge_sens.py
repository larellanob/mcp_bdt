#!/usr/bin/python

#import scipy.stats as stats
from sensitivity.systematics import *
from sensitivity.limits import get_limit

# runs 1+2+3
data_pot = 1.5e21
data_eve = 3081362.

use_systematics = False
use_sig_syst = False
if not use_systematics:
    use_sig_syst = False
use_flat_50p = True


#masses = [ 15.,20.,30.,50.,80.,100.,150.,200.,250.,300.,350.,400. ]
masses = [200.]
obs_limits = []
exp_limits = []

for m in masses:
    #directory = "root/MillichargeGamma3D_run1-2a-2b-3a/"
    #filename = "hist_BDT_scores_%imev.root"%m
    directory = "./"
    filename  = "test.root"
    '''

    get_limit() function takes as input a file with two histograms
    and one tree:

    - sig_hist and bkg_hist, histograms, correspond to the bdt scores
      of the signal and background respectively for the test/data
      sample

    - total_pot: ttree which contains two (double type) branches with
      one entry: tot_pot and tot_evt, the total pot of the simulated
      signal, and the number of simulated background events, these
      will enter into the normalization

    '''
    o_l, e_l = get_limit(directory,filename,m,use_systematics)
    obs_limits.append(o_l)
    exp_limits.append(e_l)

print("# mass, obs limits, exp limits")
for i in range(0,len(masses)):
    print(masses[i],obs_limits[i],exp_limits[i])
    
epsilon = 0.001
nhits   = 2

print("15. 1.0")
for i in range(0,len(masses)):
    sens = epsilon * pow(obs_limits[i],1./(2.+nhits*2))
    print(masses[i],sens)
print("400.0 1.0")
