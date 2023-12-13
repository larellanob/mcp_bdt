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

use_systematics = False
use_sig_syst = False
if not use_systematics:
    use_sig_syst = False
use_flat_50p = True

def get_CL_range(observations, model,unbounded_bounds,cl):

    alpha = round(1-0.01*cl,2)
    low_end  = 1000
    high_end = 2000
    obs_limit = 2000
    iterations = 0
    delta_iter = 1000
    while obs_limit >= high_end:
        iterations += 1
        if iterations > 1:
            delta_iter = iterations*10000
            high_end   += delta_iter
        poi_values = np.linspace(low_end,high_end,3)
        obs_limit, exp_limits = pyhf.infer.intervals.upper_limits.upper_limit(
            observations,model,poi_values,level=alpha,par_bounds=unbounded_bounds
        )
        print(iterations,"iterations")
    #return np.linspace(low_end,2*high_end,10)
    #lo_range = max(1000,obs_limit-obs_limit/2.)
    lo_range = 1000.
    if obs_limit > 10000:
        lo_range = obs_limit-obs_limit/2.
    up_range = obs_limit+obs_limit/2.
    return np.linspace(lo_range,up_range,50)
    

def detector_systematics(sig_or_bkg, mass,detvars="detvars"):
    # detector variations only make sense compared to CV
    # so need to do everything in relation to that
    detvars_list = [] # return this
    if detvars == "detvars":
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
    else:
        detvars_names = ["quadrature"]

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

def detector_systematics_quadrature(mass):
    # calls previous, returns both
    sig_syst = detector_systematics("signal",mass,detvars="quadrature")
    bkg_syst = detector_systematics("background",mass,detvars="quadrature")
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
    pot_tree = uprooted['total_pot']
    npot = pot_tree["tot_pot"].array()[0]
    nevt = pot_tree["tot_evt"].array()[0]
    print("POT:",npot)
    print("events:",nevt)

    pot_runs123 = 1.5e21
    eve_runs123 = 3081362.
    sig_scale, bkg_scale = pot_runs123/float(npot),eve_runs123/float(nevt)
    print(sig_scale,bkg_scale)

    # we got rid of non-pot carrying trees
    #sig_scale, bkg_scale = get_scaling(mass)

    print("\nData from the root file (before scaling): (s,b,b_err)")
    sig_bins = sig.values().tolist()
    bkg_bins = bkg.values().tolist()
    err_bins = bkg.errors().tolist()

    '''
    sig_bins = sig_bins[:25]
    bkg_bins = bkg_bins[:25]
    err_bins = err_bins[:25]
    '''
    
    print(sig_bins)
    print(bkg_bins)
    print(err_bins)

    sig_bins = [ i * sig_scale for i in sig_bins ]
    bkg_bins = [ i * bkg_scale for i in bkg_bins ]

    sig_bins = [ i * 10./3. for i in sig_bins ]
    bkg_bins = [ i * 10./3. for i in bkg_bins ]
    
    err_bins = [ math.sqrt(i) for i in bkg_bins ] # statisical error
    
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
    #sig_syst, bkg_syst = detector_systematics_sigbkg(mass) # 8 hours for same result
    if use_systematics:
        sig_syst, bkg_syst = detector_systematics_quadrature(mass)
        

        for dv in sig_syst:
            # when looping over a dict you're looping over the keys
            if dv == 'CV':
                continue
            data = model.spec['channels'][0]['samples'][0]['data']
            absolute = []
            for i in range(0,len(data)):
                absolute.append(abs(data[i]*sig_syst[dv][i]))
            new_entry = {'name':dv+"_sig",'type':'shapesys', 'data': absolute}
            #new_entry = {'name':dv+"_sig",'type':'shapesys', 'data': sig_syst[dv]}
            if use_sig_syst:
                model.spec['channels'][0]['samples'][0]['modifiers'].append(new_entry)
                print("using signal systematics")

        for dv in bkg_syst:
            print(dv)
            # when looping over a dict you're looping over the keys
            if dv == 'CV':
                continue
            data = model.spec['channels'][0]['samples'][1]['data']
            absolute = []
            for i in range(0,len(data)):
                print(data[i],bkg_syst[dv][i],data[i]*bkg_syst[dv][i])
                if use_flat_50p:
                    absolute.append(abs(data[i]*0.5))
                else:
                    absolute.append(abs(data[i]*bkg_syst[dv][i]))
                    print("BKG syst[i]:",i,bkg_syt[dv][i])
            new_entry = {'name':dv+"_bkg",'type':'shapesys', 'data': absolute}
            model.spec['channels'][0]['samples'][1]['modifiers'].append(new_entry)

    print("dumping the entire json model")
    print(json.dumps(model.spec,indent=2))

    # remodel including systematics
    if use_systematics:
        model = pyhf.Model(model.spec)
        
    #exit()

    nbins = model.config.channel_nbins['singlechannel']
    print("\nNumber of bins:",nbins)
    print(type(nbins))

    # suggested initial parameters
    init_pars = model.config.suggested_init()
    print("\nInitial parameters:",init_pars)
    print("\nSuggested bounds:",model.config.suggested_bounds())

    unbounded_bounds = model.config.suggested_bounds()
    unbounded_bounds[model.config.poi_index] = (0, 160000000)

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


    print("\nObservations")
    print(observations)

    lp = model.logpdf(pars=bkg_pars,data=observations)
    print("\nLog pdf of the model given parametrs and data")
    print(type(lp))
    print(lp)

    
    # need pyhf version >= 0.7.0 for pyhf.infer.intervals.upper_limits
    print("\npyhf version:",pyhf.__version__)

    # confidence level
    cl = 90.
    alpha = round(1-0.01*cl,2)
    #poi_values = np.linspace(0.1, 5, 50)
    #poi_values = np.linspace(5000, 10000000, 50)
    # np.linspace(start, stop, number_of_evenly_spaced_samples=50)


    #  get_CL_range()
    if use_systematics == True:
        poi_values = np.linspace(1000, 300000, 50)
        if mass == 150:
            poi_values = np.linspace(1000, 500000, 50)
        if mass == 200:
            poi_values = np.linspace(1000, 150000, 50)
        if mass == 300:
            poi_values = np.linspace(1000, 350000, 50)
        if mass == 350:
            poi_values = np.linspace(5000, 800000, 50)
        if mass == 400:
            poi_values = np.linspace(5000, 20000000, 50)
    elif use_systematics == False:
        #poi_values = np.linspace(1000, 7000,500)
        poi_values = np.linspace(1000, 12000,100)
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

    poi_values = get_CL_range(observations,model,unbounded_bounds,cl)
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
        if use_systematics:
            out_figname_pdf = "%s/limit_mass_%i_cl_%i.pdf"%(model_tag,mass,cl)
            out_figname_png = "%s/limit_mass_%i_cl_%i.png"%(model_tag,mass,cl)
        elif not use_systematics:
            out_figname_pdf = "%s/limit_mass_%i_cl_%i_no-systematics.pdf"%(model_tag,mass,cl)
            out_figname_png = "%s/limit_mass_%i_cl_%i_no-systematics.png"%(model_tag,mass,cl)
        plt.savefig(out_figname_pdf)
        plt.savefig(out_figname_png)        
    if show:
        print("\nShowing plot")
        plt.show()
    return obs_limit, exp_limits

masses = [ 15.,20.,30.,50.,80.,100.,150.,200.,250.,300.,350.,400. ]
obs_limits = []
exp_limits = []
#mass = masses[4]

#detvars_dict = detector_systematics("signal",200)
print("Printing list")
#for dv in detvars_dict:
#    print(dv,detvars_dict[dv],"\n")
#exit()

for m in masses:
    directory = "root/MillichargeGamma3D_run1-2a-2b-3a/"
    filename = "hist_BDT_scores_%imev.root"%m
    o_l, e_l = get_limit(directory,filename,m)
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
