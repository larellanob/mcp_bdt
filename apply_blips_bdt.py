#!/usr/bin/python3

import xgboost as xgb
xgb.set_config(verbosity=2)
from  matplotlib import pyplot as plt
import uproot
import pandas as pd
import os
import numpy as np
import json

# two types of file need BDT applied
## validation: for validation signal samples (overtraining check), need separated signal and bkg
### apply_scores_validation()
## not validation: for data, samples, they will not have "sig_blip_pairs" trees, etc.
### apply_scores()

def apply_scores_signal(model, sig_mass, combined_masses = False):
    # load model
    bdt = xgb.Booster()
    bdt.load_model(f"models/{model}_model.json")

    # file to validate
    input_file=f"/exp/uboone/data/users/arellano/bdt_mcp/blipreco/ntuple_241012_500kev_2hits/validation/blippairs_MillichargeBlipAna_500kev_2hits_{sig_mass}mev.root"
    if  combined_masses:
        input_file=f"/exp/uboone/data/users/arellano/bdt_mcp/blipreco/ntuple_241012_500kev_2hits/training/blippairs_MillichargeBlipAna_500kev_2hits_combinedmasses.root"

    uproot_file = uproot.open(input_file)
    sig_tree = uproot_file['sig_blip_pairs']
    bkg_tree = uproot_file['bkg_blip_pairs']
    mix_tree = uproot_file['mix_blip_pairs']
    all_branches = sig_tree.keys()
    training_vars = all_branches
    rse_vars = [ "run", "subrun", "event", "bp_blipid1", "bp_blipid2" ]
    for v in rse_vars:
        # variables not used in training
        # prediction vars must match the training
        # so we remove these, and will readd them
        # so we store them in rse dataframes
        training_vars.remove(v)
    sig_events = sig_tree.arrays(training_vars,library="pd")
    bkg_events = bkg_tree.arrays(training_vars,library="pd")
    mix_events = mix_tree.arrays(training_vars,library="pd")
    sig_rse = sig_tree.arrays(rse_vars,library="pd")
    bkg_rse = bkg_tree.arrays(rse_vars,library="pd")
    mix_rse = mix_tree.arrays(rse_vars,library="pd")

    # xgb.predict() needs DMatrix type as input
    xgb_sig = xgb.DMatrix(sig_events)
    xgb_bkg = xgb.DMatrix(bkg_events)
    xgb_mix = xgb.DMatrix(mix_events)
    samples = [(xgb_sig,'validation_sig',sig_rse),
               (xgb_bkg,'validation_bkg',bkg_rse),
               (xgb_mix,'validation_mix',mix_rse),
           ]
    if combined_masses:
        samples = [(xgb_sig,'training_sig',sig_rse),
                   (xgb_bkg,'training_bkg',bkg_rse),
                   (xgb_mix,'training_mix',mix_rse),
                   ]
        
    # output root file
    if not os.path.exists(f"root/{model}"):
        os.mkdir(f"root/{model}")
    output_root_file = f"root/{model}/validation_BDT_scores_{sig_mass}mev.root"
    if combined_masses:
        output_root_file = f"root/{model}/training_BDT_scores_combinedmasses.root"
    print("Exporting file "+output_root_file)
    with uproot.recreate(output_root_file) as f:
        for sample in samples:
            print("Obtaining predictions for sample",sample[1])
            # sample[0] is the DMatrix
            # sample[1] is the string (output tree name)
            # sample[2] is rse dataframes
            prediction = bdt.predict(sample[0],output_margin=True)
            
            # branches which will go in the file are input from dictionary
            branches = {}
            branches['bdt']         = prediction # bdt score
            branches['run']         = sample[2]['run']
            branches['subrun']      = sample[2]['subrun']
            branches['event']       = sample[2]['event']
            branches['bp_blipid1']  = sample[2]['bp_blipid1']
            branches['bp_blipid2']  = sample[2]['bp_blipid2']
            
            f[sample[1]] = branches
        print("Finished obtaining predictions")
        if not combined_masses:
            f['total_pot'] = uproot_file['total_pot'].arrays(["tot_pot"],library="pd")



gen_th = "500"
nhits  = "2"
model = "2415_blipreco_v2" #

apply_scores_validation(model,0,True) # combined_masses
exit()
masses = [ 15, 20, 30, 50, 80, 100, 150, 200, 250, 300, 350, 400 ]
for m in masses:
    apply_scores_signal(model,m)


