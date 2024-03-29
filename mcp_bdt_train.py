#!/usr/bin/python3

import xgboost as xgb
xgb.set_config(verbosity=2)
from  matplotlib import pyplot as plt
import uproot
import pandas
import os
import numpy as np
import numpy.ma as ma # masked numpy arrays, to remove invalid data
from sklearn.preprocessing import StandardScaler
import json

debug = False
necessary_dirs = ["models","root","img"]
for d in necessary_dirs:
    if not os.path.exists(d):
        os.mkdir(d)



fcl_th = "def-th"
gen_th = "5000"
nhits  = "1"

#tag   = "%s_%skev_%shits"%(fcl_th,gen_th,nhits)
#tag = "MillichargeGamma3D_run1-2a-2b-3a"
tag = 'MillichargeBlip_%skev'%(gen_th)
# model = "%s_230424"%(tag)
model = "%s_231106"%(tag)
input_root_file = 'root/%s/preselection_%skev_%shits_combinedmass.root'%(tag,gen_th,nhits)
#################################################################################
#################################################################################
#################################################################################
######### UPROOT 
# open root file, grab tree
try:
    uproot_file = uproot.open(input_root_file)#["tsig_train"]
except OSError:
    print("Can't find input root file. Terminating")
    print(input_root_file)
    exit()


if debug:
    print("\nLet's have a look inside the uprooted file\n")
    print(uproot_file)
    print(uproot_file.keys())

    print("\nLook inside a tree within the file\n")
    print(uproot_file['tsig_train'])
    print(uproot_file['tsig_train'].keys())
    print(uproot_file['tsig_train'].keys()[0])

    print("\nLook inside a branch within a tree within the file\n")
    print(uproot_file['tsig_train']['px1'])
    print(uproot_file['tsig_train'].arrays()) # needs awkward for this

    print("\nThis is the official way to look into branches\n")

tree = uproot_file['tsig_train']
branches = tree.arrays()


if debug:
    print("branches",branches)
    print(branches['px1'])
    print("\nThis is the official way to look into events\n")
    print(branches['px1'][0:5])
    print("\nLooks good?\n")

    print("####################################################################")
    print("                              featmap                               ")
    print("####################################################################")

# make featmap file
featmap_vars = {
        'int':[  'nhits1', 'nhitsO1', 'nhitsM1',  'nhits2',
            'nhitsO2', 'nhitsM2', 'dca_nh2_behind', 'dca_nh2_mid', 'dca_nh2_ahead', 'dca_nhO_behind',
            'dca_nhO_mid', 'dca_nhO_ahead', 'dca_nhM_behind', 'dca_nhM_mid',
                 'dca_nhM_ahead'],
#            'sum_nhits_1','sum_nhits_2','sum_nhits'],
        'q':[ 'px1', 'py1', 'pz1', 'log10_adc1', 'px2', 'py2', 'pz2', 'log10_adc2',
            'dist', 'theta', 'phi', 'sqrt_dca_behind',
            'sqrt_dca_ahead' ],
        'i':[],
        }

# each row a feature (tree branch)
# columns: counter (integer), feature name, type (q,int)
featmap_counter = 0

featmap_file_name = "models/%s_featmap.txt"%(model)
featmap_file = open(featmap_file_name,'w')
featmap_list = []
for i in tree.keys():
    for j in featmap_vars:
        if i in featmap_vars[j]:
            line = str(featmap_counter)+"\t"+str(i)+"\t"+str(j)+"\n"
            featmap_file.write(line)
            if debug:
                print("--------")
                print(i,j)
                print(line)
                print("--------")
            featmap_list.append(i)
            featmap_counter+=1
featmap_file.close()
# do not close! we will keep using this
# yes close, then we open again

if debug:
    print(featmap_list)
    print(len(featmap_list))
    print("####################################################################")
    print("                               train                                ")
    print("####################################################################")


train_sig = uproot_file['tsig_train'].arrays(featmap_list,library="np")
train_bkg = uproot_file['tbg_train'].arrays(featmap_list,library="np")
test_sig = uproot_file['tsig_test'].arrays(featmap_list,library="np")
test_bkg = uproot_file['tbg_test'].arrays(featmap_list,library="np")



#train_sig_concat = np.stack((train_sig['px1'],train_sig['px2']))
train_sig_data = np.stack(list(train_sig[y] for y in featmap_list))
train_bkg_data = np.stack(list(train_bkg[y] for y in featmap_list))
test_sig_data  = np.stack(list(test_sig[y]  for y in featmap_list))
test_bkg_data  = np.stack(list(test_bkg[y]  for y in featmap_list))

print("TRAIN SIG DATA MAX MIN",np.amax(train_sig_data),np.amin(train_sig_data))

# full train dataset
train_data = np.concatenate((train_sig_data,train_bkg_data),1)
# remove numbers which are too large (and unphysical)
train_data = ma.masked_greater(train_data,1e90) # mask numbers above 1e90
train_data = train_data.filled(-10000) # actually replace them with -10000
# do the same for other datasets

train_sig_data = ma.masked_greater(train_sig_data,1e90)
train_sig_data = train_sig_data.filled(-10000)
train_bkg_data = ma.masked_greater(train_bkg_data,1e90)
train_bkg_data = train_bkg_data.filled(-10000)
test_sig_data  = ma.masked_greater(test_sig_data,1e90)
test_sig_data  = test_sig_data.filled(-10000)
test_bkg_data  = ma.masked_greater(test_bkg_data,1e90)
test_bkg_data  = test_bkg_data.filled(-10000)



if debug:
    print("train_sig_concat")
    print(train_sig_data)
    print("train_sig_data dim, shape")
    print(train_sig_data.ndim,train_sig_data.shape)
    print("list of features and first entry")
    for y in featmap_list:
        print(y,train_sig[y][0])


    print("Pandas version:",pandas.__version__)
    print("xgb version:",xgb.__version__)

    print('len train sig',len(train_sig),'bkg',len(train_bkg))
    print('len test sig',len(test_sig),'bkg',len(test_bkg))
    print("\nAll good so far\n")


#truth_label = np.concatenate((np.zeros(len(train_bkg_data[0])),np.ones(len(train_sig_data[0]))))
truth_label = np.concatenate((np.ones(len(train_sig_data[0])),np.zeros(len(train_bkg_data[0]))))
if debug:
    print(truth_label)
    print(truth_label.size)
    print(type(truth_label))
truth_label=truth_label.ravel()
if debug:
    print("flattened")
    print(truth_label)
    print(truth_label.size)
    print(type(truth_label))
train_data = np.transpose(train_data)
if debug:
    print("train_data dim, shape")
    print(train_data.ndim,train_data.shape)

    print("\nsig_data\n")
    print(train_sig_data.size)
    print(train_sig_data.ndim,train_sig_data.shape)
l =[1]*len(train_sig_data[0])
l2=[0]*len(train_bkg_data[0])
l3=[1]*len(test_sig_data[0])
l4=[0]*len(test_bkg_data[0])

print("Len of l (signal train) and first five elements::")
print(len(l))
print(l[:5])

train_sig_data = np.transpose(train_sig_data)
train_bkg_data = np.transpose(train_bkg_data)
test_sig_data = np.transpose(test_sig_data)
test_bkg_data = np.transpose(test_bkg_data)

### training proper
xgb_train     = xgb.DMatrix(train_data,label=truth_label)
xgb_train_sig = xgb.DMatrix(train_sig_data)
xgb_train_sig.set_label(l)

xgb_train_bkg = xgb.DMatrix(train_bkg_data)
xgb_train_bkg.set_label(l2)
xgb_test_sig  = xgb.DMatrix(test_sig_data)
xgb_test_sig.set_label(l3)
xgb_test_bkg  = xgb.DMatrix(test_bkg_data)
xgb_test_bkg.set_label(l4)

watchlist = [(xgb_train_sig,'train_sig'), # you are training with this
             (xgb_train_bkg,'train_bkg'), # so watching metrics not useful?
             (xgb_test_sig,'test_sig'),
             (xgb_test_bkg,'test_bkg')]
if debug:
    print("xgb dmatrices made, watchlist made successfully")

# https://xgboost.readthedocs.io/en/stable/python/python_api.html#module-xgboost.training
#xgboost.train(params,               #Booster params.
#              dtrain,               #(DMatrix) – Data to be trained.
#              num_boost_round=10,   #(int) – Number of boosting iterations.

#              evals=(),             # (list of pairs (DMatrix, string))
                                     # List of validation sets for which metrics
                                     # will evaluated during training.
                                     # Validation metrics will help us track
                                     # the performance of the model.
#              obj=None,
#              feval=None,
#              maximize=None,
#              early_stopping_rounds=None,
#              evals_result=None,
#              verbose_eval=True,
#              xgb_model=None,
#              callbacks=None
#              )

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
#    ('eval_metric','auc'), # i added this
#    ('eval_metric','map'), # i added this
    ('eval_metric','error'), # i added this
    ('learning_rate',0.1)
]

# put a high number here and "early_stopping_rounds" will cut the
# training short if there's no improvement in X rounds
num_rounds = 300
#num_rounds = 1

evals_result={}
print("Training begin")
bdt = xgb.train(parameters,xgb_train,num_rounds,evals=watchlist,evals_result=evals_result,early_stopping_rounds=20)
#bdt = xgb.train(parameters,xgb_train,num_rounds)
print("Training complete")
print(bdt.best_iteration)
#exit()

print(json.dumps(evals_result,indent=2))
for i in evals_result:
    #print(i)
    #print(evals_result[i]["logloss"])

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    logloss1 = ax1.plot(evals_result[i]["logloss"],'b-')
    error1 = ax2.plot(evals_result[i]["error"],'r-')
    plt.title(i)
    ax1.set_xlabel("iteration")
    ax1.set_ylabel("logloss",color='b')
    ax2.set_ylabel("error",color='r')
    plt.show()


#exit()
# pawel saves multiple runs, look into that
bdt.save_model("models/%s_model.json"%(model))
bdt.dump_model("models/%s_dump.txt"%(model))
# I don't know why it's dumping with feature names already


if debug:
    print("####################################################################")
    print("                              results                               ")
    print("####################################################################")

#xgb.plot_importance(bdt)
#plt.show()

#xgb.plot_tree(bdt)
#plt.show()

featmap_file = open(featmap_file_name,'r')
if debug:
    print("scores:")
    print("featmap file type",type(featmap_file_name))
    print("featmap file",featmap_file_name)

scores = bdt.get_score(fmap=featmap_file_name,importance_type='total_gain')
print("keys:")
keys = sorted(scores.keys(),key = lambda k: scores[k])

print("for keys:")
for k in keys:
    print(k,scores[k])


if debug:
    print("####################################################################")
    print("                               export                               ")
    print("####################################################################")


# we now want to make a graph for checking overtraining
# currently saving a root file with combined masses, then using a root macro to plot
# want to just plot directly here and don't bother creating that root file
# the root files we do need are mass by mass, not combined
    
output_root_file = "root/%s/%s_BDT_scores.root"%(tag,model)
print("Exporting file "+output_root_file)


samples = watchlist
print("wlist")
samples[0] += (uproot_file['tsig_train'].arrays(["run","evt"],library="pd"),)
print("passed 0")
samples[1] += (uproot_file['tbg_train'].arrays(["run","evt"],library="pd"),)
print("passed 1")
samples[2] += (uproot_file['tsig_test'].arrays(["run","evt"],library="pd"),)
print("passed 2")
samples[3] += (uproot_file['tbg_test'].arrays(["run","evt"],library="pd"),)
print("defined samples")
with uproot.recreate(output_root_file) as f:
    print("with file as f")
    for sample in samples:
        print("sample")
        # sample[0] is the DMatrix
        # sample[1] is the string
        #prediction = bdt.predict(sample[0], iteration_range=(0,num_rounds),output_margin=True)
        prediction = bdt.predict(sample[0],output_margin=True)
        
        # branches which will go in the file are input from dictionary
        # they get added backwards for some reason
        branches = {}
        branches['run'] = sample[2]['run'].to_numpy()
        branches['evt'] = sample[2]['evt'].to_numpy()
        branches['bdt'] = prediction

        f[sample[1]] = branches
print("masses")
masses = [15,20,30,50,80,100,150,200,250,300,350,400]
out_test_files = []
for m in masses:
    print("MASS:",m)

    #output_test_root_file = "root/%s/%s_BDT_scores_mass%s.root"%(tag,model,m)
    #test_file = 'root/%s/preselection_%skev_%shits_%smev.root'%(tag,gen_th,nhits,m)
    test_file = 'root/%s/preselection_%smev.root'%(tag,m)
    print("opened file",test_file)
    test_uproot = uproot.open(test_file)

    print(test_uproot)
    
    bdt_test = xgb.Booster()
    bdt_test.load_model("models/%s_model.json"%(model))

    sig_test = test_uproot['tsig_test'].arrays(featmap_list,library='np')
    bkg_test = test_uproot['tbg_test'].arrays(featmap_list,library='np')

    sig_test = np.stack(list(sig_test[y] for y in featmap_list))
    bkg_test = np.stack(list(bkg_test[y] for y in featmap_list))

    sig_test = np.transpose(sig_test)
    bkg_test = np.transpose(bkg_test)

    sig_test = ma.masked_greater(sig_test,1e90)
    sig_test = sig_test.filled(-10000)
    bkg_test = ma.masked_greater(bkg_test,1e90)
    bkg_test = bkg_test.filled(-10000)

    xgb_test_sig2 = xgb.DMatrix(sig_test)
    xgb_test_bkg2 = xgb.DMatrix(bkg_test)


    samples = [ (xgb_test_sig2,'test_sig'),
                (xgb_test_bkg2,'test_bkg')]
    # tuple adding (appending as last element)
    samples[0] += (test_uproot['tsig_test'].arrays(["run","evt"],library="np"),)
    samples[1] += (test_uproot['tbg_test'].arrays(["run","evt"],library="np"),)


    out_test_file = 'root/%s/test_BDT_scores_%smev.root'%(tag,m)
    out_test_files.append(out_test_file)
    print("saving file",out_test_file)
    with uproot.recreate(out_test_file) as f2:
        for sample in samples:
            # sample[0] is the DMatrix
            # sample[1] is the string
            prediction = bdt.predict(sample[0], output_margin=True)

            # sample[2] are the new elements, run and evt
            # not a part of the training but might be useful when testing
            # branches which will go in the file are input from dictionary
            # they get added backwards for some reason
            branches = {}
            branches['run'] = sample[2]['run']
            branches['evt'] = sample[2]['evt']
            branches['bdt'] = prediction
            
            f2[sample[1]] = branches
        f2['total_pot'] = test_uproot['total_pot'].arrays(["tot_pot","tot_evt"],library="np")

            
print(out_test_files)
# Make_BDT_histograms.cxx already loops over masses
for tf in out_test_files:
    root_cmd = 'root -l -b -q macro/Make_BDT_histograms.cxx\'("%s")\''%(tf)
    os.system(root_cmd)


# make overtraining check figures
model_number = model.split('_')[-1]
root_cmd = 'root -l -b -q macro/overtraining_check.cxx\'("%s","%s")\''%(tag,model_number)
print("running")
print(root_cmd)
os.system(root_cmd)
