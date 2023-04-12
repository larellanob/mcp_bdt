#!/usr/bin/python3

import xgboost as xgb
xgb.set_config(verbosity=2)
from  matplotlib import pyplot as plt
import uproot
import pandas
import os
import numpy as np
import numpy.ma as ma # masked numpy arrays, to remove invalid data

fcl_th = "def-th"
gen_th = "1000"
nhits  = "2"
mass   = "350"
tag_mass  = "%s_%skev_%shits_%smev"%(fcl_th,gen_th,nhits,mass)
tag   = "%s_%skev_%shits"%(fcl_th,gen_th,nhits)
model = "%s_230301"%(tag)
detvar_dir = "systematics/detvar_%s_v08_00_00_61"%(tag_mass)

detvars = [ 'CV',
            'wiremod_x',
            'wiremod_yz',
            'wiremod_anglexz',
            'wiremod_angleyz',
            'wiremod_dEdx',
            'ly_down',
            'ly_rayleigh',
            'ly_atten',
            'sce',
            'recomb'
           ]

featmap_file_name = "models/%s_featmap.txt"%(model)
featmap_file = open(featmap_file_name,'r')
featmap_list = []
for line in featmap_file:
    split = line.split()
    featmap_list.append(split[1])

print(featmap_list)

for detvar in detvars:
    detvar_file = '%s/preselection_%s.root'%(detvar_dir,detvar)
    detvar_uproot = uproot.open(detvar_file)
    bdt_test = xgb.Booster()
    bdt_test.load_model("models/%s_model.json"%(model))

    sig_test = detvar_uproot['tsig_test'].arrays(featmap_list,library='np')
    bkg_test = detvar_uproot['tbg_test'].arrays(featmap_list,library='np')

    # correct dimensions
    sig_test = np.stack(list(sig_test[y] for y in featmap_list))
    bkg_test = np.stack(list(bkg_test[y] for y in featmap_list))
    sig_test = np.transpose(sig_test)
    bkg_test = np.transpose(bkg_test)

    # mask to remove infs and nans
    sig_test = ma.masked_greater(sig_test,1e90)
    sig_test = sig_test.filled(-10000)
    bkg_test = ma.masked_greater(bkg_test,1e90)
    bkg_test = bkg_test.filled(-10000)

    
    xgb_test_sig2 = xgb.DMatrix(sig_test)
    xgb_test_bkg2 = xgb.DMatrix(bkg_test)

    #l_bkg=[0]*len(bkg_trans)
    #l_sig=[1]*len(sig_trans)
    #xgb_test_bkg2.set_label(l_bkg)
    #xgb_test_sig2.set_label(l_sig)

    watchlist = [
        (xgb_test_sig2,'test_sig',detvar_uproot['tsig_test'].arrays(["run","evt"],library="pd")),
        (xgb_test_bkg2,'test_bkg',detvar_uproot['tbg_test'].arrays(["run","evt"],library="pd"))
    ]
                   
    out_test_file = '%s/test_BDT_scores_%s.root'%(detvar_dir,detvar)
    print("saving file",out_test_file)

    with uproot.recreate(out_test_file) as f2:
        for sample in watchlist:
            # sample[0] is the DMatrix
            # sample[1] is the string
            prediction = bdt_test.predict(sample[0], output_margin=True)
            
            # branches which will go in the file are input from dictionary
            # they get added backwards for some reason
            branches = {}
            branches['run'] = sample[2]['run'].to_numpy()
            branches['evt'] = sample[2]['evt'].to_numpy()
            branches['bdt'] = prediction
            
            f2[sample[1]] = branches

    # make bdt root histograms
    os.system('root -l -b -q macro/Make_BDT_histograms.cxx\'("%s","%s")\''%(out_test_file,detvar))

print("Getting detvar ratios")
os.system('root -l -b -q macro/Get_detvar_ratios.cxx\'("%s")\''%detvar_dir)
