import matplotlib.pyplot as plt
import numpy as np
import pyhf
from pyhf.contrib.viz import brazil
import math
import pathlib
import json
import uproot
from math import sqrt
import time

def get_CL_range(mass,observations, model,unbounded_bounds,cl,fast_mode = False):
    '''This function obtains a quick range of values for the parameter
    of interest (POI) mu by performing upper_limits.upper_limit on
    just three points, bound by low_end and high_end, such that
    desired value of mu "obs_limit" is contained within these bounds
    but not to close to either border. What usually happens is that
    the desired mu is much larger than high_end and thus upper_limit
    defaults to put its value in the upper bound, so we will
    continuously raise high_end until obs_limit is lower than it

    '''
    debug = True
    if debug:
        print("______________________ debugging get_CL_range -______________________")
        print(mass)
        print(observations)
        print(model)
        print(unbounded_bounds)
        print(cl)
    alpha = round(1-0.01*cl,2)
    low_end  = 0.00001 # lower bound of mu (POI)
    high_end = 0.01 # upper bound of mu (POI)
    obs_limit = 20 # value of mu at which CLs(mu) = 0.1
    '''
    if mass < 70:
        high_end = 0.1
    elif mass >= 70 and mass < 150:
        high_end = 1.0
    elif mass >= 150 and mass < 300:
        high_end = 20.0
        obs_limit = 30.0
    elif mass >= 300 and mass < 400:
        high_end = 200.0
        obs_limit = 250.0
    '''
    iterations = 0
    delta_iter = 0.01
    if debug:
        print("about to enter while")
        print("Initial: high end =",high_end,"obs limit =",obs_limit)
    while obs_limit >= high_end:
        iterations += 1
        if iterations > 3:
            delta_iter = delta_iter*1.5
            high_end   += delta_iter
        elif iterations > 10:
            delta_iter = iterations*10
            high_end   += delta_iter
        else:
            high_end += delta_iter
        poi_values = np.linspace(low_end,high_end,3)
        #print("about to call up lim")
        obs_limit, exp_limits = pyhf.infer.intervals.upper_limits.upper_limit(
            observations,model,poi_values,level=alpha,par_bounds=unbounded_bounds
        )
        if debug:
            print(iterations,"iterations. low end =",low_end,", high end =",high_end,"obs limit =",obs_limit,"delta iter =",delta_iter)
    if debug:
        print("passed while")
        print("Final: high end =",high_end,"obs limit =",obs_limit)
    #return np.linspace(low_end,2*high_end,10)
    #lo_range = max(1000,obs_limit-obs_limit/2.)
    lo_range = low_end
    if obs_limit > 10000:
        lo_range = obs_limit-obs_limit/2.
    up_range = obs_limit+obs_limit/2.
    if not fast_mode:
        return np.linspace(lo_range,up_range,50)
    else:
        return np.linspace(lo_range,up_range,5)


def get_limit(directory,filename,mass,timestamp,errors_mode="stat-only",verbose = False,fast_mode = False):
    if timestamp == "":
        print("ERROR: No valid timestamp. Exiting.")
        exit()
    print(timestamp)
    print(f"\nSensitivity for mass {mass:.1f}")
    

    #input_root_file = '../root/post_bdtcut_s_b_mass_%i.root'%(mass)
    input_root_file = directory+filename
    uprooted = uproot.open(input_root_file)

    h_sig  =  uprooted["h_logit_20_sig"   ].values().tolist()
    h_elas =  uprooted["h_logit_01_elas"  ].values().tolist()
    h_kdar =  uprooted["h_logit_02_kdar"  ].values().tolist()
    h_bkg  =  uprooted["h_logit_03_bg"    ].values().tolist()
    h_cos  =  uprooted["h_logit_04_cosmic"].values().tolist()
    h_ext  =  uprooted["h_logit_05_ext"   ].values().tolist()

    # figure out after which point historgams contain only zeroes
    last_nonzero_bin_allhists = 0
    for h in [ h_sig, h_elas, h_kdar, h_bkg, h_cos, h_ext]:
        last_nonzero_bin = 0
        for b in range(0,len(h)):
            if h[b] != 0.0 and last_nonzero_bin < b:
                last_nonzero_bin += 1
            last_nonzero_bin_allhists = max(last_nonzero_bin_allhists,last_nonzero_bin)

    h_sig  = h_sig[ :last_nonzero_bin_allhists+1]
    h_elas = h_elas[:last_nonzero_bin_allhists+1]
    h_kdar = h_kdar[:last_nonzero_bin_allhists+1]
    h_bkg  = h_bkg[ :last_nonzero_bin_allhists+1]
    h_cos  = h_cos[ :last_nonzero_bin_allhists+1]
    h_ext  = h_ext[ :last_nonzero_bin_allhists+1]

    # make sure especially that signal doesn't have non-zero bins
    for i in range(0,len(h_sig)):
        if h_sig[i] == 0.0:
            print("Found zero bin in signal!")
            h_sig[i] = 1.0e-5
    #h_sig[8]=0.00015594538000000016
            
        
    err_sig  = uprooted["h_logit_20_sig"   ].errors().tolist()[:last_nonzero_bin_allhists+1]
    err_elas = uprooted["h_logit_01_elas"  ].errors().tolist()[:last_nonzero_bin_allhists+1]
    err_kdar = uprooted["h_logit_02_kdar"  ].errors().tolist()[:last_nonzero_bin_allhists+1]
    err_bkg  = uprooted["h_logit_03_bg"    ].errors().tolist()[:last_nonzero_bin_allhists+1]
    err_cos  = uprooted["h_logit_04_cosmic"].errors().tolist()[:last_nonzero_bin_allhists+1]
    err_ext  = uprooted["h_logit_05_ext"   ].errors().tolist()[:last_nonzero_bin_allhists+1]
    # stat error sqrt(n)
    if verbose == True:
        print("Signal/bkg values histograms/lists")
        for i in [ h_sig, h_elas, h_kdar, h_bkg, h_cos, h_ext]:
            print(len(i))
            print(i)
        print("Signal/bkg errors histograms/lists")
        for err in [ err_sig, err_elas, err_kdar, err_bkg, err_cos, err_ext]:
            print("err",err)
            print(err)

    l_all_bkg =     [ h_elas,h_kdar,h_bkg,h_cos,h_ext ]
    l_all_bkg_err = [ err_elas,err_kdar,err_bkg,err_cos,err_ext ]

    model_dict = {
        "channels": [
            {
                "name": "singlechannel",
                "samples":  [
                    {
                        "name": "signal",
                        "data": h_sig,
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
                        "data": h_cos,
                        "modifiers": [
                            {
                                "name": "genie_err",
                                "type": "staterror", # all of these used to be staterror
                                "data": err_cos
                            }
                        ]
                    },
                    {
                        "name": "nubkg",
                        "data": h_bkg,
                        "modifiers": [
                            {
                                "name": "genie_err",
                                "type": "staterror",
                                "data": err_bkg
                            }
                        ]
                    },
                    {
                        "name": "kdar",
                        "data": h_kdar,
                        "modifiers": [
                            {
                                "name": "genie_err",
                                "type": "staterror",
                                "data": err_kdar
                            }
                        ]
                    },
                    {
                        "name": "ext",
                        "data": h_ext,
                        "modifiers": [
                            {
                                "name": "ext_err",
                                "type": "staterror",
                                "data": err_ext
                            }
                        ]
                    },
                    {
                        "name": "elas",
                        "data": h_elas,
                        "modifiers": [
                            {
                                "name": "elas_err",
                                "type": "staterror",
                                "data": err_elas
                            }
                        ]
                    }
                ] #end channel sample
            } 
        ] # end channels
    } # end model dict
    #exit()

    model = pyhf.Model(model_dict)
    if verbose == True:
        print("successfully modeled, dumping model")
        print(model)
        print(json.dumps(model.spec, indent=2))

    # ERRORS MODES
    # stat-only (skips all of these)
    # shapesys (detector variations using quadrature and shapesys)
    # histosys (detector variations using quadrature and histosys (not implemented yet?))
    if errors_mode == "shapesys":
        syst_file = "/gluster/data/microboone/millicharge/wirecell/detvars/400mev_1-15_quadrature.root"
        syst_uprt = uproot.open(syst_file)
        quad = syst_uprt["h_quad"].values().tolist()

        # we need to convert this "ratio" uncertainty into an absolute one
        data = model.spec['channels'][0]['samples'][0]['data'] #sample 0 is signal
        absolute = []
        for i in range(0,len(data)):
            absolute.append(abs(data[i]*quad[i]*200.0))
        # give it appropriate dictionary format to add to the model
        new_entry = {'name':'quadrature_30mev','type':'shapesys', 'data': absolute}
        # add systematics as modifier
        model.spec['channels'][0]['samples'][0]['modifiers'].append(new_entry)
        if verbose:
            print("using Systematics as shapesys")
            print("quad:")
            print(quad)
            print("added systematics in quadrature, before and after:")
            print(data)
            print(absolute)

    elif errors_mode == "histosys":
        syst_file = "/gluster/data/microboone/millicharge/wirecell/detvars/400mev_1-15_quadrature.root"
        syst_uprt = uproot.open(syst_file)
        quad = syst_uprt["h_quad"].values().tolist()

        # we need to convert this "ratio" uncertainty into an absolute one
        data = model.spec['channels'][0]['samples'][0]['data'] #sample 0 is signal
        hi_data = []
        lo_data = []
        for i in range(0,len(data)):
            hi_data.append(data[i]+data[i]*quad[i])
            lo_data.append(data[i]-data[i]*quad[i])
            #lo_data.append(max(0.0,data[i]-data[i]*quad[i]))
        # give it appropriate dictionary format to add to the model
        new_entry = {'name':'quadrature_30mev_histosys','type':'histosys', 'data': {"hi_data": hi_data, "lo_data": lo_data}}
        # add systematics as modifier
        model.spec['channels'][0]['samples'][0]['modifiers'].append(new_entry)
        if verbose:
            print("using HISTOSYS")
            print("quad:")
            print(quad)

            print("                 ======== MODEL EXPECTED DATA= ================")
            print(model.expected_data)
            print(f"  aux data: {model.config.auxdata}")
        #print(f"      down: {model.expected_data([-1.0])}")
        #print(f"   nominal: {model.expected_data([0.0])}")
        #print(f"        up: {model.expected_data([1.0])}")
        
    # suggested initial parameters
    init_pars = model.config.suggested_init()
    # extend bounds
    unbounded_bounds = model.config.suggested_bounds()
    unbounded_bounds[model.config.poi_index] = (0, 1000000000000000000.0) # set large number for weak signals
    if verbose == True:
        print("\nSuggested bounds:",model.config.suggested_bounds())
        print("\nInitial parameters:",init_pars)
        print("\nAuxdata")
        print(model.config.auxdata)
        print(model)

    # one histogram bkg_bins combining sum of all bkgs
    bkg_bins = []
    for b in zip(*l_all_bkg):
        if verbose == True:
            print("backgrounds values entries and length")
            print(b)
            print(len(b))
        sum_of_bkgs=0
        for bi in range(0,len(b)):
            sum_of_bkgs+=b[bi]
        bkg_bins.append(sum_of_bkgs)

    #same thing but for errors, and save that to model auxdata
    bkg_errbins = []
    for b in zip(*l_all_bkg_err):
        if verbose == True:
            print("background error entries and length")
            print(b)
            print(len(b))
        sum_of_bkg_err=0
        for bi in range(0,len(b)):
            sum_of_bkg_err+=b[bi]
        bkg_errbins.append(sum_of_bkg_err)

    print(json.dumps(model.spec, indent=2)) 
        
    if verbose == True:
        print("model config auxdata",model.config.auxdata)
        print("bkg_errbins",bkg_errbins)
        print(model.config.auxdata)
        print(len(model.config.auxdata))

    if verbose == True:
        num = 0.0
        denom = 0.0
        for i,j in zip(l_all_bkg,l_all_bkg_err):
            num+= j[0]
            denom+= i[0]
        num = sqrt(num)
        print(num/denom)
        print(model.config.parameters)
        print(model.config.param_set)
        print("model config auxdata\n",model.config.auxdata)


    # bkg only hypothesis
    observations = bkg_bins + model.config.auxdata

    if verbose == True:
        print(bkg_bins)
        print(model.config.auxdata)
        print(observations)
        # need pyhf version >= 0.7.0 for pyhf.infer.intervals.upper_limits
        print("\npyhf version:",pyhf.__version__)


    cl = 90.
    alpha = round(1-0.01*cl,2)
    if not fast_mode:
        poi_values = get_CL_range(mass,observations,model,unbounded_bounds,cl)
    else:
        poi_values = get_CL_range(mass,observations,model,unbounded_bounds,cl,fast_mode=True)
    
    if verbose == True:
        print("##############################################################")
        print("######### Calling upper_limits:")
        print("##############################################################")
        print(f"\nh_sig (len {len(h_sig)})")
        print(h_sig)
        print(f"\nobservations (len {len(observations)})")
        print(observations)
        print(f"\npoi values (len {len(poi_values)})")
        print(poi_values)
        print(f"\nunbounded_bounds (len {len(unbounded_bounds)})")
        print(unbounded_bounds)
        print(f"\nmodel:")
        print(model)
        print(json.dumps(model.spec, indent=2)) 

    
    obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(
        observations, model, poi_values, level=alpha, return_results=True, par_bounds=unbounded_bounds
    )

    fig, ax = plt.subplots()
    fig.set_size_inches(10.5, 7)
    ax.set_title(f"Hypothesis Tests, mass {mass}, μ ={obs_limit:.4f}")

    artists = brazil.plot_results(poi_values, results, test_size=alpha, ax=ax)


    save = True
    show = False
    #model_tag = "230123_normal_thresholds"
    model_tag = directory.replace("root","img")

    if save:
        print(f"\nSaving plot to {model_tag}/{timestamp}")
        pathlib.Path(model_tag+"/"+timestamp).mkdir(parents=True,exist_ok=True)
        out_figname_pdf = f"{model_tag}/{timestamp}/limit_mass_{mass}_cl_{cl}_{errors_mode}.pdf"
        out_figname_png = f"{model_tag}/{timestamp}/limit_mass_{mass}_cl_{cl}_{errors_mode}.png"
        plt.savefig(out_figname_pdf)
        plt.savefig(out_figname_png)        
    if show:
        print("\nShowing plot")
        plt.show()
    return obs_limit, exp_limits

all_masses = [15.,20.,30.,50.,80.,100.,150.,200.,250.,300.,350.,400.]
meson  = 'rho'
meson = ""
masses = all_masses
if meson == 'pi0':
    masses = [15.,20., 30., 50.]
if meson == 'eta':
    masses = [15.,20., 30., 50., 80., 100., 150., 200., 250.]
if meson == 'etp':
    #masses = all_masses
    #masses = [300., 350.,400.]
    #masses = [250.]
    #masses = [400.]
    '''
    i'm having lots of problems with etp mass 400mev

    the upper_limit function is saying it's dividing by zero. this
    converged for the previous model (combined masses, non-separate
    per parents) and the only difference is the signal histogram
    bins. they are not even that different.

    the function converges by manually setting one magic bin to the
    one from the previous model/sample: see commented line
    h_sig[8]=0.00015594538000000016 above
    
    '''
    masses = [15.,20., 30., 50., 80., 100., 150., 200., 250., 300., 350.]
if meson == 'rho':
    masses = [15.,20., 30., 50., 80., 100., 150., 200., 250., 300., 350.]
    
obs_limits = []
exp_limits = []

model = '20000kev_1hits_combinedmasses'
#model = '20000kev_combinedmass_per_parent_'+meson

fastmode = False # for testing purposes
masses = all_masses
#masses = [15.,20.,30.,50.]
masses = [80.,100.,150.]
# errors_mode that do something:
# stat-only - no systematic uncertainties
# shapesys
# histosys
errors_mode = "histosys"
if errors_mode != "stat-only" and errors_mode != "shapesys" and errors_mode != "histosys":
    print("ERROR: errors mode not valid")

t = time.localtime()
timestamp = time.strftime('d%Y%m%d_t%H%M%S',t)

for m in masses:
    directory = f"/gluster/data/microboone/millicharge/sensitivity/{model}/"
    filename  = f"{model}_{int(m)}mev.root"
    
    o_l, e_l = get_limit(directory,filename,m,timestamp,fast_mode=fastmode,errors_mode=errors_mode)
    obs_limits.append(o_l)
    exp_limits.append(e_l)

print("# mass, obs limits, exp limits")
for i in range(0,len(masses)):
    print(masses[i],obs_limits[i],exp_limits[i])
    
epsilon = 0.001
nhits   = 1

sens_txt_file_out = f"sensitivity/mysensitivities/sens_{timestamp}_{errors_mode}.txt"
output_file = open(sens_txt_file_out,"w")
output_file.write("{} 1.0\n".format(masses[0]))
for i in range(0,len(masses)):
    sens = epsilon * pow(obs_limits[i],1./(2.+nhits*2))
    output_file.write("{} {}\n".format(masses[i],sens))
output_file.write("{} 1.0\n".format(masses[-1]))
