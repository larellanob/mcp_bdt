#!/usr/bin/python

import json
import matplotlib.pyplot as plt
import numpy as np
import pyhf
from pyhf.contrib.viz import brazil
import uproot
import scipy.stats as stats
import math
import pathlib

# runs 1+2+3
data_pot = 1.5e21
data_eve = 3081362.


# normal thresholds
def get_scaling(mass):
    # returns sig scaling, bkg scaling
    eve = 0
    pot = 0
    pot_runs123 = 1.5e21
    eve_runs123 = 3081362.
    if mass == 100:
        eve = 11617
        pot = 1.48e27
    elif mass == 150:
        eve = 14335
        pot = 2.44e27
    elif mass == 200:
        eve = 11310
        pot = 2.56e27
    elif mass == 300:
        eve = 9902
        pot = 4.76e27
    elif mass == 350:
        eve = 9950
        pot = 9.76e27
    elif mass == 400:
        eve = 9953
        pot = 2.28e29
    else:
        return 0
    return pot_runs123/float(pot),eve_runs123/float(eve)

'''
# low thresholds

def get_scaling(mass):
    # returns sig scaling, bkg scaling
    eve = 0
    pot = 0
    pot_runs123 = 1.5e21
    eve_runs123 = 3081362.
    if mass == 100:
        eve = 12012
        pot = 1.52e27
    elif mass == 150:
        eve = 9550
        pot = 1.63e27
    elif mass == 200:
        eve = 9600
        pot = 2.16e27
    elif mass == 300:
        eve = 9500
        pot = 4.70e27
    elif mass == 350:
        eve = 9553
        pot = 9.36e27
    elif mass == 400:
        eve = 9950
        pot = 2.15e29
    else:
        return 0
    return pot_runs123/float(pot),eve_runs123/float(eve)
'''

def detector_systematics(sig_or_bkg, mass):
    # detector variations only make sense compared to CV
    # so need to do everything in relation to that
    detvars_list = [] # return this
    detvars_names = [
        "CV",
        "wiremod_x",
        "wiremod_yz",
        "wiremod_anglexz",
        "wiremod_angleyz",
        "wiremod_dEdx",
        "ly_down",
        "ly_rayleigh",
        "ly_atten",
        "sce",
        "recomb"
        ]

    detvars_dict = {}
    
    ratios_root_file = 'systematics/detvar_def-th_1000kev_2hits_%imev_v08_00_00_61/Quadrature_%s.root'%(mass,sig_or_bkg)

    ur = uproot.open(ratios_root_file)
    for dv in detvars_names:
        try:
            detvars_list.append(ur[dv].values().tolist())
            detvars_dict[dv] = ur[dv].values().tolist()
        except:
            print("No histogram found for detvar",dv)
            detvars_dict[dv] = []

    # returns dictionary of samples : histograms (list of bin values)
    return detvars_dict

def detector_systematics_sigbkg(mass):
    # calls previous, returns both
    sig_syst = detector_systematics("signal",mass)
    bkg_syst = detector_systematics("background",mass)
    return sig_syst,bkg_syst

def get_limit(directory, filename, mass):
    print("\nSensitivity for mass %.1f"%(mass))

    #input_root_file = '../root/post_bdtcut_s_b_mass_%i.root'%(mass)
    input_root_file = directory+filename
    uprooted = uproot.open(input_root_file)

    print(uprooted)
    print(uprooted.keys())

    sig = uprooted['sig_hist']
    bkg = uprooted['bkg_hist']

    # from nupy arrays to python list
    # would be better if you made pyhf read numpy arrays?
    # pyhf.set_backend("numpy")

    #sig_scale = get_sig_scaling(mass)
    #bkg_scale = get_bkg_scaling(mass)
    sig_scale, bkg_scale = get_scaling(mass)

    print("\nData from the root file (before scaling): (s,b,b_err)")
    sig_bins = sig.values().tolist()
    bkg_bins = bkg.values().tolist()
    err_bins = bkg.errors().tolist()

    print(sig_bins)
    print(bkg_bins)
    print(err_bins)

    sig_bins = [ i * sig_scale for i in sig_bins ]
    bkg_bins = [ i * bkg_scale for i in bkg_bins ]

    sig_bins = [ i * 10./3. for i in sig_bins ]
    bkg_bins = [ i * 10./3. for i in bkg_bins ]
    
    err_bins = [ math.sqrt(i) for i in bkg_bins ]
    

    '''
    use only first n bins
    n = 5
    sig_bins = sig_bins[:n]
    bkg_bins = bkg_bins[:n]
    err_bins = err_bins[:n]
    '''
    # empty bins are being a problem
    '''
    print("\nSETTING ALL EMPTY BKG BINS TO 1")
    for i in range(len(bkg_bins)):
        if bkg_bins[i] == 0.:
            bkg_bins[i] = 1.
            err_bins[i] = 1.
   '''

    print("\nData from the root file: (s,b,b_err)")
    print(sig_bins)
    print(bkg_bins)
    print(err_bins)


    model = pyhf.simplemodels.uncorrelated_background(
        signal=sig_bins,
        bkg=bkg_bins,
        bkg_uncertainty=err_bins
    )
    
    print("\nWe have a model of type")
    print(type(model))

    print("dumping the entire json model")
    print(json.dumps(model.spec,indent=2))

    ## systematics
    sig_syst, bkg_syst = detector_systematics_sigbkg(mass)

    for dv in sig_syst:
        # when looping over a dict you're looping over the keys
        if dv == 'CV':
            continue
        data = model.spec['channels'][0]['samples'][0]['data']
        absolute = []
        for i in range(0,len(data)):
            absolute.append(data[i]*sig_syst[dv][i])
        #new_entry = {'name':dv,'type':'shapesys', 'data': absolute}
        new_entry = {'name':dv,'type':'shapesys', 'data': sig_syst[dv]}
        model.spec['channels'][0]['samples'][0]['modifiers'].append(new_entry)

    for dv in bkg_syst:
        # when looping over a dict you're looping over the keys
        if dv == 'CV':
            continue
        data = model.spec['channels'][0]['samples'][1]['data']
        absolute = []
        for i in range(0,len(data)):
            absolute.append(data[i]*sig_syst[dv][i])
        #new_entry = {'name':dv,'type':'shapesys', 'data': absolute}
        new_entry = {'name':dv,'type':'shapesys', 'data': sig_syst[dv]}
        model.spec['channels'][0]['samples'][1]['modifiers'].append(new_entry)
    
    print("dumping the entire json model")
    print(json.dumps(model.spec,indent=2))

    #exit()

    nbins = model.config.channel_nbins['singlechannel']
    print("\nNumber of bins:",nbins)
    print(type(nbins))

    # suggested initial parameters
    init_pars = model.config.suggested_init()
    print("\nInitial parameters:",init_pars)
    print("\nSuggested bounds:",model.config.suggested_bounds())

    unbounded_bounds = model.config.suggested_bounds()
    unbounded_bounds[model.config.poi_index] = (0, 20000000)

    # background only expectation!
    # sig=0, bkg=1
    #model.expected_actualdata([0.0]+nbins*[1.0])
    bkg_pars = init_pars.copy()
    # signal (parameter of interest) set to 0
    bkg_pars[model.config.poi_index] = 0
    print("\nExpected actual data (bkg only)")
    print(model.expected_actualdata(bkg_pars))

    print("\nAuxdata")
    print(model.config.auxdata)

    # observations = background only
    observations = bkg_bins + model.config.auxdata

    '''
    # randomize observations?
    observations = []
    for i in range(0,len(bkg_bins)):
        pois = np.random.poisson(bkg_bins[i])
        print(bkg_bins[i],pois)
        observations.append(pois)
    print(observations)

    observations += model.config.auxdata
    '''
    
    print("\nObservations")
    print(observations)

    lp = model.logpdf(pars=bkg_pars,data=observations)
    print("\nLog pdf of the model given parametrs and data")
    print(type(lp))
    print(lp)

    print("\nFit")
    #np.seterr(divide = 'ignore')
    fit = pyhf.infer.mle.fit(data=observations,pdf=model)
    print(fit)
    print("\u03BC\u0302 = ",fit[0])
    print("\u03B3\u0302 = ",fit[1:])
    # should I scale mu (signal strength) at this point?

    # CLs
    print("\nRunning CLs")
    CLs_obs, CLs_exp = pyhf.infer.hypotest(
        1.0, # null hypothesis "BSM physics process exists" (????)
        # in pyhf documentation they call null hypothesis when you have full signal (mu = 1)
        observations, 
        model,
        test_stat = "qtilde", # can use q or qtilde (?)
        return_expected_set=True,
    )

    print(model.config.suggested_bounds()[model.config.poi_index])

    print(f"      Observed CLs: {CLs_obs:.4f}")
    for expected_value, n_sigma in zip(CLs_exp, np.arange(-2, 3)):
        print(f"Expected CLs({n_sigma:2d} σ): {expected_value:.4f}")

    # need pyhf version >= 0.7.0 for pyhf.infer.intervals.upper_limits
    print("\npyhf version:",pyhf.__version__)
    
    # confidence level
    cl = 90.
    alpha = round(1-0.01*cl,2)
    #poi_values = np.linspace(0.1, 5, 50)
    #poi_values = np.linspace(5000, 10000000, 50)

    poi_values = np.linspace(1000, 100000, 100)
    if mass == 200:
        poi_values = np.linspace(1000, 150000, 100)
    if mass == 300:
        poi_values = np.linspace(1000, 350000, 100)
    if mass == 350:
        poi_values = np.linspace(5000, 800000, 100)
    if mass == 400:
        poi_values = np.linspace(5000, 20000000, 100)

    # comment for testing
    #return
    
    print("\nGetting limit for mass",mass)
    obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(
        observations, model, poi_values, level=alpha, return_results=True, par_bounds=unbounded_bounds
    )
    print(f"Upper limit (obs): μ = {obs_limit:.4f}")
    print(f"Upper limit (exp): μ = {exp_limits[2]:.4f}")


    fig, ax = plt.subplots()
    fig.set_size_inches(10.5, 7)
    ax.set_title("Hypothesis Tests, mass %i, μ =%.4f"%(mass,obs_limit))

    artists = brazil.plot_results(poi_values, results, test_size=alpha, ax=ax)

    save = True
    show = False
    #model_tag = "230123_normal_thresholds"
    model_tag = directory.replace("root","img")
    if save:
        print("\nSaving plot")
        pathlib.Path(model_tag).mkdir(parents=True,exist_ok=True)
        plt.savefig("%s/limit_mass_%i_cl_%i.pdf"%(model_tag,mass,cl))
    if show:
        print("\nShowing plot")
        plt.show()
    return obs_limit, exp_limits

masses = [ 100.,150.,200.,300.,350.,400. ]
obs_limits = []
exp_limits = []
#mass = masses[4]

detvars_dict = detector_systematics("signal",200)
print("Printing list")
for dv in detvars_dict:
    print(dv,detvars_dict[dv],"\n")
#exit()

for m in masses:
    directory = "root/def-th_1000kev_2hits/"
    filename = "hist_BDT_scores_1000kev_2hits_%imev.root"%m
    o_l, e_l = get_limit(directory,filename,m)
    obs_limits.append(o_l)
    exp_limits.append(e_l)

print("# mass, obs limits, exp limits")
for i in range(0,len(masses)):
    print(masses[i],obs_limits[i],exp_limits[i])
    
epsilon = 0.001
nhits   = 2

print("100.0 1.0")
for i in range(0,len(masses)):
    sens = epsilon * pow(obs_limits[i],1./(2.+nhits*2))
    print(masses[i],sens)
print("400.0 1.0")
