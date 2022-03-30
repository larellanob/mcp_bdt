#!/usr/bin/python3

import xgboost as xgb
from  matplotlib import pyplot as plt
import uproot
import pandas
import os

debug = True
necessary_dirs = ["models","root","img"]
for d in necessary_dirs:
    if not os.path.exists(d):
        os.mkdir(d)

model = "my_model_2"
input_root_file = 'root/spacepoints.root_output.root'

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
        'q':[ 'px1', 'py1', 'pz1', 'adc1', 'px2', 'py2', 'pz2', 'adc2',
            'dist', 'theta', 'phi', 'dca_behind', 'dca_mid',
            'dca_ahead', 'dca_adc_behind', 'dca_adc_mid', 'dca_adc_ahead',
            'mindist_closest', 'mindist_adc_closest' ],
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
            featmap_list.append(i)
            featmap_counter+=1
#featmap_file.close()
# do not close! we will keep using this

# i have created a featmap.txt file
# is this all pawel's make_featmap do?
# do i need to store this somewhere else? (a list, etc)

if debug:
    print("####################################################################")
    print("                               train                                ")
    print("####################################################################")

# equivalent to pawel's train.py

# he gets using params = params_getter(args.model)
# this function looks into the file but i'm not sure what returns

train_sig = uproot_file['tsig_train'].arrays(featmap_list,library="pd")
train_bkg = uproot_file['tbg_train'].arrays(featmap_list,library="pd")
test_sig = uproot_file['tsig_test'].arrays(featmap_list,library="pd")
test_bkg = uproot_file['tbg_test'].arrays(featmap_list,library="pd")



if debug:
    print(train_bkg)
    print("train_bkg type:",type(train_bkg))
    print(featmap_list)
    reduced = train_bkg[0:5]
    print(reduced)
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


print("Pandas version:",pandas.__version__)
print("xgb version:",xgb.__version__)


if debug:
    print('len train sig',len(train_sig),'bg',len(train_bkg))
    print('len test sig',len(test_sig),'bg',len(test_bkg))
    print("\nAll good so far\n")

truth_label = [1]*len(train_sig) + [0]*len(train_bkg)

if debug:
    print(truth_label)
    print(pandas.concat([train_sig,train_bkg]))


### training proper

# why these four?
xgb_train     = xgb.DMatrix(pandas.concat([train_sig,train_bkg]),truth_label)
xgb_train_sig = xgb.DMatrix(train_sig,label=[1]*len(train_sig))
xgb_train_bkg = xgb.DMatrix(train_bkg,label=[0]*len(train_bkg))
xgb_test_sig  = xgb.DMatrix(test_sig,label=[1]*len(test_sig))
xgb_test_bkg  = xgb.DMatrix(test_bkg,label=[0]*len(test_bkg))

watchlist = [(xgb_train_sig,'train_sig'), # you are training with this
             (xgb_train_bkg,'train_bkg'), # so watching metrics not useful?
             (xgb_test_sig,'test_sig'),
             (xgb_test_bkg,'test_bkg')]

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

xgb.plot_importance(bdt)
#plt.show()

xgb.plot_tree(bdt)
#plt.show()

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
