#!/usr/bin/python3

import xgboost as xgb
xgb.set_config(verbosity=2)
from  matplotlib import pyplot as plt
import uproot
import pandas as pd
import os
import numpy as np
import json
from sklearn.model_selection import train_test_split



gen_th = "500"
nhits  = "2"
model = "2415_blipreco_v2" #

#input_root_file = "/exp/uboone/data/users/arellano/bdt_mcp/blipreco/ntuple_241012_500kev_2hits/training/blippairs_MillichargeBlipAna_500kev_2hits_combinedmasses.root" # v1
input_root_file = "/exp/uboone/data/users/arellano/bdt_mcp/blipreco/numi/v1/run1/blippairs/sig_training/hist_blippairs_numi_run1_sig_combinedmasses.root" #v2
######### UPROOT 
# open root file, grab tree
try:
    uproot_file = uproot.open(input_root_file)
except OSError:
    print("Can't find input root file. Terminating")
    print(input_root_file)
    exit()

"""
file contains 3 trees: 
- sig_blip_pairs Both blips signal
- bkg_blip_pairs Both blips background
- mix_blip_pairs One blip signal, one blip background

what we'll do

1. define variables used for training
   - get them from the files themselves
   - exclude things not to be used (event, run, subrun, etc.)
   - you need featmap_list list for
     train_sig = uproot_file['tsig_train'].arrays(featmap_list,library="np")

2. define parameters
   - optimize

3. train using sig and bkg only, not mix
   - use 50/50 train/test
   - use many metrics and save the plots (logloss/auc)
   - use early stopping rounds

"""

# we open the tree and get its list of branches (variables)
# branches from sig only, bkg must have the same always (?)
sig_tree = uproot_file['sig_blip_pairs']
bkg_tree = uproot_file['bkg_blip_pairs']
all_branches = sig_tree.keys()

# we train with all variables except dont_train_list
train_vars = all_branches
dont_train_list = [ "run", "subrun", "event", "bp_blipid1", "bp_blipid2" ]
for v in dont_train_list:
    train_vars.remove(v)

# split events used for training/testing
training_frac = 0.8
testing_frac = 1.0-training_frac

sig_events = sig_tree.arrays(train_vars,library="pd")
#sig_events = sig_events[0:10] # test with fewer events
sig_train, sig_test = train_test_split(sig_events,
                                       test_size=testing_frac)
bkg_events = bkg_tree.arrays(train_vars,library="pd")
#bkg_events = bkg_events[0:10] # test with fewer events
bkg_train, bkg_test = train_test_split(bkg_events,
                                       test_size=testing_frac)
full_train = pd.concat([sig_train,bkg_train])
full_test = pd.concat([sig_test,bkg_test])

# make xgboost dmatrices
# separate for metrics, then one combined
xgb_train = xgb.DMatrix(full_train)
xgb_sig_train = xgb.DMatrix(sig_train)
xgb_bkg_train = xgb.DMatrix(bkg_train)
xgb_test = xgb.DMatrix(full_test)
xgb_sig_test = xgb.DMatrix(sig_test)
xgb_bkg_test = xgb.DMatrix(bkg_test)
# set "truth" labels
xgb_train.set_label([1]*xgb_sig_train.num_row()+[0]*xgb_bkg_train.num_row())
xgb_sig_train.set_label([1]*xgb_sig_train.num_row())
xgb_bkg_train.set_label([0]*xgb_bkg_train.num_row())
xgb_test.set_label([1]*xgb_sig_test.num_row()+[0]*xgb_bkg_test.num_row())
xgb_sig_test.set_label([1]*xgb_sig_test.num_row())
xgb_bkg_test.set_label([0]*xgb_bkg_test.num_row())

print("Number sig train",xgb_sig_train.num_row())
print("Number bkg train",xgb_bkg_train.num_row())
print("Number sig test",xgb_sig_test.num_row())
print("Number bkg test",xgb_bkg_test.num_row())

# training parameters
parameters = [
    ('booster', 'dart'),
    ('rate_drop', 0.1),
    ('skip_drop', 0.5),
    ('max_depth', 6),
    ('eta', 0.3),
    ('objective', 'binary:logistic'),
    ('subsample', 0.5),
#    ('tree_method', 'hist'),
    ('tree_method', 'auto'),
    ('scale_pos_weight', None),
    ('eval_metric','logloss'), # i added this
    ('eval_metric','auc'), # i added this
    ('eval_metric','error'), # i added this
    ('learning_rate',0.1)
]
num_rounds = 1000
num_early_stop = 20

# samples to evaluate
watchlist = [(xgb_train,'Train'),
             (xgb_test,'Test')]

# empty dict where eval metrics will be stored
evals_result = {}
# TRAINING

print("Training begin")
bdt = xgb.train(parameters,
                xgb_train,
                num_rounds,
                evals=watchlist,
                evals_result=evals_result,
                early_stopping_rounds=num_early_stop
                )
print("Finished training!")

# save the model and dump it
bdt.save_model(f"models/{model}_model.json")
bdt.dump_model(f"models/{model}_dump.txt")

# print the dumped dump the model
#print(json.dumps(evals_result,indent=2))

# plot the eval metrics
'''
for i in evals_result:
    #print(i)
    #print(evals_result[i]["logloss"])

    os.makedirs(f"img/{model}", exist_ok=True)

    metrics = ["logloss","error","auc"]
    colours = ["blue","red","green"]
    mtitles = ["Log Loss","Error","AUC"]
    
    # we will (try to) do all in the same figure
    # so we need many axes
    ax_counter = 0
    for m in zip(metrics,colours,mtitles):
        ax_counter+=1
        met = m[0]
        col = m[1]
        tit = m[2]

        # single 
        fig, ax = plt.subplots()
        plt.title(i)
        pl = ax.plot(evals_result[i][met],col)
        ax.set_xlabel("Iteration")
        ax.set_ylabel(tit,color=col)
        ax.figure.savefig(f"img/{model}/{i}_{met}.pdf")
        
        # multiple thing
        if ax_counter == 1:
            # if we're on the first axis, make the base figure
            multifig, multiax = plt.subplots()
            plt.title(i)
            multiax.plot(evals_result[i][met],col)
            multiax.set_xlabel("Iteration")
            multiax.set_ylabel(tit,color=col)
            multiax.figure.savefig(f"img/{model}/{i}.pdf")
        else:
            # else clone the first axis
            multiax1 = multiax.twinx()
            multiax1.plot(evals_result[i][met],col)
            if ax_counter > 2:
                # if we're above two axis, keep moving it to the right by 50
                multiax1.spines['right'].set_position(('outward',50*(ax_counter-2)))
            multiax1.set_ylabel(tit,color=col)
            # need bbox_inches='tight' to "recalculate" the width
            # otherwise last axes fall out of the pdf
            multiax1.figure.savefig(f"img/{model}/{i}.pdf",bbox_inches='tight')
'''
metrics = ["logloss","error","auc"]
mtitles = ["Log Loss","Error","AUC"]

for m in zip(metrics,mtitles):
    plt.clf()
    plt.plot(evals_result['Train'][m[0]],color='b',label="Train")
    plt.plot(evals_result['Test'][m[0]], color='r',label="Test")
    plt.title("XGBoost training - Blip pairs")
    plt.xlabel("Iteration")
    plt.ylabel(m[1])
    plt.legend()
    print("saving file to",f"img/{model}/{m[0]}.pdf")
    os.makedirs(f"img/{model}", exist_ok=True)
    plt.savefig(f"img/{model}/{m[0]}.pdf")

                
# output root file
if not os.path.exists(f"root/{model}"):
    os.mkdir(f"root/{model}")
output_root_file = f"root/{model}/training_BDT_scores.root"
print("Exporting file "+output_root_file)

# samples to export
samples = [(xgb_sig_train,'train_sig'),
           (xgb_bkg_train,'train_bkg'),
           (xgb_sig_test,'test_sig'),
           (xgb_bkg_test,'test_bkg')
           ]

# export
with uproot.recreate(output_root_file) as f:
    for sample in samples:
        # sample[0] is the DMatrix
        # sample[1] is the string
        prediction = bdt.predict(sample[0],output_margin=True)

        # branches which will go in the file are input from dictionary
        # they get added backwards for some reason
        branches = {}
        # you can add more branches e.g.
        #branches['bp_px1'] = sample[1]["bp_px1"] # not tested!
        branches['bdt'] = prediction

        f[sample[1]] = branches
        
