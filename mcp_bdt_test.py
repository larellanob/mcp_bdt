#!/usr/bin/python3
import sys
import xgboost as xgb
import pkg_resources
print("Using xgboost version",pkg_resources.get_distribution("xgboost").version)
from  matplotlib import pyplot as plt
import uproot
import pandas
import numpy as np
import numpy.ma as ma # masked numpy arrays, to remove invalid data
from sklearn.preprocessing import StandardScaler

#test_sample = "neutrino_overlay_270_files"
#test_sample = "1000kev_2hits_400mev"
#train_model="model_221207_combinedmass_pot"
train_model="model_230123_normal_thresholds"
#test_sample = "1nu1cos"
#test_sample = "neutrino"
#test_sample = "100_normal_thresholds"

print("length of arguments:",len(sys.argv))
if len(sys.argv) == 1:
    print("no argument, using default sample:",test_sample)
if len(sys.argv) == 2:
    test_sample = "1000kev_2hits_%smev"%(sys.argv[1])
    print(test_sample)
elif len(sys.argv) == 3:
    test_sample = "1000kev_2hits_%smev_%s"%(sys.argv[1],sys.argv[2])
    print(test_sample)

#input_root_file = 'root/spacepoints.root_output_invertedADC.root'
input_root_file = "root/%s/preselection_%s.root"%(train_model,test_sample)
uproot_file = uproot.open(input_root_file)


# check that all features in the featmap file are also in the test
# sample root file
tree = uproot_file['tsig_test']
fmap_sample = tree.keys() # list of features
# featmap file is associated to the training model
featmap_file_name = "models/%s_featmap.txt"%(train_model)
featmap_list = []
with open(featmap_file_name,'r') as f:
    lines = f.readlines()
featmap_counter = 0

for l in range(0,len(lines)):
    sp = lines[l].split('\t')
    featmap_list.append(sp[1])
    featmap_counter += 1

for fmap_model in featmap_list:
    if fmap_model not in fmap_sample:
        print("ERROR:",fmap_model,"not in test sample!")
        exit()
    
'''
featmap_vars = {
    'int':['nhits1', 'nhitsO1', 'nhitsM1',  'nhits2',
           'nhitsO2', 'nhitsM2', 'dca_nh2_behind', 'dca_nh2_mid', 'dca_nh2_ahead', 'dca_nhO_behind',
           'dca_nhO_mid', 'dca_nhO_ahead', 'dca_nhM_behind', 'dca_nhM_mid',
           'dca_nhM_ahead',
           'sum_nhits_1','sum_nhits_2','sum_nhits'],
    'q':[ 'px1', 'py1', 'pz1', 'adc1', 'px2', 'py2', 'pz2', 'adc2',
          'dist', 'theta', 'phi', 'dca_behind', 'dca_mid',
          'dca_ahead', 'dca_adc_behind', 'dca_adc_mid', 'dca_adc_ahead',
          'mindist_closest', 'mindist_adc_closest' ],
    'i':[],
}


featmap_diff = []
temp_counter = featmap_counter
for i in tree.keys():
    for j in featmap_vars:
        if i in featmap_vars[j]:
            line = str(featmap_counter)+"\t"+str(i)+"\t"+str(j)+"\n"
            featmap_file.write(line)
            featmap_list.append(i)
            featmap_counter+=1
    if temp_counter == featmap_counter:
        featmap_diff.append(i)
    temp_counter = featmap_counter
'''

bdt = xgb.Booster()
bdt.load_model("models/%s_model.json"%(train_model))


#t_test_sig = uproot_file['tsig_test'].arrays(featmap_list,library="pd",filter_branch=lambda b: b.dca.mid < 100000)
t_test_sig = uproot_file['tsig_test'].arrays(featmap_list,library="np")
t_test_bkg = uproot_file['tbg_test'].arrays(featmap_list,library="np")
test_sig = np.stack((t_test_sig[y] for y in featmap_list))
test_bkg = np.stack((t_test_bkg[y] for y in featmap_list))
test_sig = np.transpose(test_sig)
test_sig = ma.masked_greater(test_sig,1e90)
test_sig = test_sig.filled(-10000)
test_bkg = np.transpose(test_bkg)
test_bkg = ma.masked_greater(test_bkg,1e90)
test_bkg = test_bkg.filled(-10000)
xgb_test_bkg  = xgb.DMatrix(test_bkg)
xgb_test_sig  = xgb.DMatrix(test_sig)
l_bkg=[0]*len(test_bkg)
l_sig=[1]*len(test_sig)
xgb_test_bkg.set_label(l_bkg)
xgb_test_sig.set_label(l_sig)

# scaled sig and bkg (remove nan and inf)
#scaler = StandardScaler()
#t_test_bkg_scaled = scaler.fit_transform(t_test_bkg)
#t_test_sig_scaled = scaler.fit_transform(t_test_sig)
#xgb_test_sig  = xgb.DMatrix(t_test_sig_scaled,label=[1]*len(t_test_sig_scaled))
#xgb_test_bkg  = xgb.DMatrix(t_test_bkg_scaled,label=[0]*len(t_test_bkg_scaled))





#########################################
#########################################
#### file output
#########################################
#########################################

output_root_file = "root/%s/test_%s_BDT_scores.root"%(train_model,test_sample)

sig_runevt = uproot_file['tsig_test'].arrays(["run","evt"],library="pd")
bkg_runevt = uproot_file['tbg_test'].arrays(["run","evt"],library="pd")

samples = [(xgb_test_sig,'test_sig',t_test_sig, sig_runevt),
           (xgb_test_bkg,'test_bkg',t_test_bkg, bkg_runevt)]

with uproot.recreate(output_root_file) as f:
    for sample in samples:
        # sample[0] is the DMatrix
        # sample[1] is the string
        prediction = bdt.predict(sample[0], output_margin=True)

        # branches which will go in the file are input from dictionary
        # they get added backwards for some reason
        branches = {}
        branches['run'] = sample[3]['run'].to_numpy()
        branches['evt'] = sample[3]['evt'].to_numpy()
        branches['bdt'] = prediction

        f[sample[1]] = branches


