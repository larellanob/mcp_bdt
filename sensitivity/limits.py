import matplotlib.pyplot as plt
import numpy as np
import pyhf
from pyhf.contrib.viz import brazil
import math
import pathlib
import json
import uproot

from sensitivity.samples import Sample

def get_CL_range(observations, model,unbounded_bounds,cl):

    alpha = round(1-0.01*cl,2)
    low_end  = 10
    high_end = 2000
    obs_limit = 2000
    iterations = 0
    delta_iter = 100
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
    
def get_limit(directory, filename, mass, use_systematics=False):
    print("\nSensitivity for mass %.1f"%(mass))

    #input_root_file = '../root/post_bdtcut_s_b_mass_%i.root'%(mass)
    input_root_file = directory+filename
    uprooted = uproot.open(input_root_file)

    print(uprooted)
    print(uprooted.keys())

    '''
    sample_list = []

    for sampletype in [ "sig", "bkg" ]:
        counter = 0
        for s in range(0,10):
            try:
                current_sample=Sample(uprooted,f'{sampletype}{counter}')
                sample_list.append(current_sample)
                #sig_list.append(uprooted[f'{sampletype}{counter}_hist'])
            except:
                print(f'found {counter} {sampletype} histograms')
                break
            else:
                print(f'\tfound hist: tree {sampletype}{counter}_hist')
                counter += 1


    # we got rid of non-pot carrying trees
    #sig_scale, bkg_scale = get_scaling(mass)

    print("\nData from the root file (before scaling): (s,b,b_err)")
    s_list = []
    s_err_list = []
    for s in sample_list:
        s_bins = s.get_hist_values()
        s_err_bins = s.get_err_values()
        s_list.append(s_bins)
        s_err_list.append(s_err_bins)
        print(s_bins)
        
        print("\nData from the root file: (s,b,b_err)")
        print(s_bins)
        #print(bkg_bins)
        #print(err_bins)
        

    print("this is what's going into simplemodels.uncorrelated_background:")
    print(s_list[0])
    print(s_list[1])
    print(s_err_list[1])
    model = pyhf.simplemodels.uncorrelated_background(
        #signal=sig_bins,
        #bkg=bkg_bins,
        #bkg_uncertainty=err_bins
        signal=s_list[0],
        bkg=s_list[1],
        bkg_uncertainty=s_err_list[1]
    )
    print(json.dumps(model.spec,indent=2))
    '''

    sig0 = Sample(uprooted,"sig0")
    bkg0 = Sample(uprooted,"bkg0") #
    bkg1 = Sample(uprooted,"bkg1") # 
    print("exiting")



    # build model "manually"

    model_dict = {
        "channels": [
            {
                "name": "singlechannel",
                "samples":  [
                    {
                        "name": "signal",
                        "data": sig0.get_hist_values(),
                        "modifiers": [
                            {
                                "name": "mu",
                                "type": "normfactor",
                                "data": None
                            }
                        ]
                    },
                    {
                        "name": "overlay",
                        "data": bkg0.get_hist_values(),
                        "modifiers": [
                            {
                                "name": "stat",
                                "type": "staterror",
                                "data": bkg0.get_err_values()
                            }
                        ]
                    },
                    {
                        "name": "nubkg",
                        "data": bkg1.get_hist_values(),
                        "modifiers": [
                            {
                                "name": "stat",
                                "type": "staterror",
                                "data": bkg1.get_err_values()
                            }
                        ]
                    }
                ] #end channel sample
            } 
        ] # end channels
    } # end model dict

    model = pyhf.Model(model_dict)
    
    #exit()
    
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
    print("\nUnbounded bounds:",unbounded_bounds)

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
    bkg_bins = [ i+j for i,j in zip(bkg0.get_hist_values(),bkg1.get_hist_values())] # one histogram combining sum of all bkgs
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

    print("getting POI values")
    #poi_values = get_CL_range(observations,model,unbounded_bounds,cl)
    poi_values = np.linspace(0.1,10,50)
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
