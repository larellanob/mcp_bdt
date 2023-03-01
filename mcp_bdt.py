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


debug = False
necessary_dirs = ["models","root","img"]
for d in necessary_dirs:
    if not os.path.exists(d):
        os.mkdir(d)

model = "model_230123_normal_thresholds"
input_root_file = 'root/%s/preselection_1000kev_2hits_combinedmass.root'%(model)
#################################################################################
#################################################################################
#################################################################################
######### UPROOT 
# open root file, grab tree
try:
    uproot_file = uproot.open(input_root_file)#["tsig_train"]
except OSError:
    print("Can't find input root file. Terminating")
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
print("branches",branches)


if debug:
    print(branches['px1'])
    print("\nThis is the official way to look into events\n")
    print(branches['px1'][0:5])
    print("\nLooks good?\n")

    print("####################################################################")
    print("                              featmap                               ")
    print("####################################################################")

# make featmap file
params = {'px1'}
featmap_vars = {
        'int':[  'nhits1', 'nhitsO1', 'nhitsM1',  'nhits2',
            'nhitsO2', 'nhitsM2', 'dca_nh2_behind', 'dca_nh2_mid', 'dca_nh2_ahead', 'dca_nhO_behind',
            'dca_nhO_mid', 'dca_nhO_ahead', 'dca_nhM_behind', 'dca_nhM_mid',
                 'dca_nhM_ahead'],
#            'sum_nhits_1','sum_nhits_2','sum_nhits'],
        'q':[ 'px1', 'py1', 'pz1', 'log10_adc1', 'px2', 'py2', 'pz2', 'log10_adc2',
            'dist', 'theta', 'phi', 'sqrt_dca_behind', 'sqrt_dca_mid',
            'sqrt_dca_ahead', 'dca_log10_adc_behind', 'dca_log10_adc_mid', 'dca_log10_adc_ahead',
            'sqrt_mindist_closest', 'mindist_log10_adc_closest' ],
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
            print("--------")
            print(i,j)
            print(line)
            print("--------")
            featmap_list.append(i)
            featmap_counter+=1
featmap_file.close()
# do not close! we will keep using this
# yes close, then we open again
print(featmap_list)
print(len(featmap_list))

if debug:
    print("####################################################################")
    print("                               train                                ")
    print("####################################################################")

# equivalent to pawel's train.py

# he gets using params = params_getter(args.model)
# this function looks into the file but i'm not sure what returns

train_sig = uproot_file['tsig_train'].arrays(featmap_list,library="np")
train_bkg = uproot_file['tbg_train'].arrays(featmap_list,library="np")
test_sig = uproot_file['tsig_test'].arrays(featmap_list,library="np")
test_bkg = uproot_file['tbg_test'].arrays(featmap_list,library="np")



if debug:
    print('train background')
    print(train_bkg)
    print("train_bkg type:")
    print(type(train_bkg))
    print("train_bkg type:")
    print(type(train_bkg))
    print("train background ['px1']")
    print(train_bkg['px1'])
    print("len train bkg ['px1']")
    print(len(train_bkg['px1']))
    print('featmap list')
    print(featmap_list)
    #reduced = train_bkg[0:5]
    #print(reduced)
    print("len featmap list",len(featmap_list))
    print("len t train sig",len(train_sig))
    
# take 70% and use for train, remaining for test
#train_fraction = 0.7
#train_sig = t_train_sig[:int(len(t_train_sig)*train_fraction)]
#test_sig = t_train_sig[int(len(t_train_sig)*train_fraction):]
#train_bkg = t_train_bkg[:int(len(t_train_bkg)*train_fraction)]
#test_bkg = t_train_bkg[int(len(t_train_bkg)*train_fraction):]

# remove nans
#train_sig = pandas.DataFrame(train_sig).fillna(0)
#test_sig  = pandas.DataFrame(test_sig).fillna(0)
#train_bkg = pandas.DataFrame(train_bkg).fillna(0)
#test_bkg  = pandas.DataFrame(test_sig).fillna(0)

#numpy train_sig

#train_sig_concat = np.stack((train_sig['px1'],train_sig['px2']))
train_sig_data = np.stack((train_sig[y] for y in featmap_list))
train_bkg_data = np.stack((train_bkg[y] for y in featmap_list))
test_sig_data = np.stack((test_sig[y] for y in featmap_list))
test_bkg_data = np.stack((test_bkg[y] for y in featmap_list))

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
test_sig_data = ma.masked_greater(test_sig_data,1e90)
test_sig_data = test_sig_data.filled(-10000)
test_bkg_data = ma.masked_greater(test_bkg_data,1e90)
test_bkg_data = test_bkg_data.filled(-10000)



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

#truth_label = [1]*len(train_sig) + [0]*len(train_bkg)


truth_label = np.concatenate((np.zeros(len(train_bkg_data[0])),np.ones(len(train_sig_data[0]))))
print(truth_label)
print(truth_label.size)
print(type(truth_label))
truth_label=truth_label.ravel()
print("flattened")
print(truth_label)
print(truth_label.size)
print(type(truth_label))
train_data = np.transpose(train_data)
print("train_data dim, shape")
print(train_data.ndim,train_data.shape)

print("\nsig_data\n")
print(train_sig_data.size)
print(train_sig_data.ndim,train_sig_data.shape)
l = [1]*len(train_sig_data[0])
l2=[0]*len(train_bkg_data[0])
l3=[1]*len(test_sig_data[0])
l4=[0]*len(test_bkg_data[0])

print(len(l))

train_sig_data = np.transpose(train_sig_data)
train_bkg_data = np.transpose(train_bkg_data)
test_sig_data = np.transpose(test_sig_data)
test_bkg_data = np.transpose(test_bkg_data)

### training proper
xgb_train     = xgb.DMatrix(train_data,label=truth_label)
xgb_train_sig = xgb.DMatrix(train_sig_data)
xgb_train_sig.set_label(l)
#xgb_train_sig = xgb.DMatrix(train_sig_data,label=[1]*len(train_sig_data[0]))
#xgb_train_bkg = xgb.DMatrix(train_bkg_data,label=[0]*len(train_bkg_data[0]))
#xgb_test_sig  = xgb.DMatrix(test_sig_data,label=[1]*len(test_sig_data[0]))
#xgb_test_bkg  = xgb.DMatrix(test_bkg_data,label=[0]*len(test_bkg_data[0]))
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
    ('max_depth', 6),
    ('eta', 0.3),
    ('objective', 'binary:logistic'),
    ('subsample', 0.5),
    ('tree_method', 'hist'),
    ('scale_pos_weight', None),
    ('rate_drop', 0.1),
    ('eval_metric','logloss'), # i added this
    ('eval_metric','error'), # i added this
    ('skip_drop', 0.5)
]

num_rounds = 30

evals_result={}
print("Training begin")
bdt = xgb.train(parameters,xgb_train,num_rounds,evals=watchlist,evals_result=evals_result)
print("Training complete")

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

output_root_file = "root/%s_BDT_scores.root"%(model)


#watchlist = [(xgb_train_sig,'train_sig'), # you are training with this
#             (xgb_train_bkg,'train_bkg'), # so watching metrics not useful?
#             (xgb_test_sig,'test_sig'),
#             (xgb_test_bkg,'test_bkg')]

# pawel saves multiple runs, look into that

print(xgb_train_sig)
print("type of train_sig is ",type(train_sig))
print(train_sig)
print(train_sig['px1'])


with uproot.recreate(output_root_file) as f:
    for sample in watchlist:
        # sample[0] is the DMatrix
        # sample[1] is the string
        prediction = bdt.predict(sample[0], output_margin=True)
        print(type(prediction))
        print(type(sample[0]))
        print(prediction.size)
        f[sample[1]] = {"bdt": prediction}
        #f[sample[1]] = {"test": 2.3}
