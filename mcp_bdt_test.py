#!/usr/bin/python3

import xgboost as xgb
from  matplotlib import pyplot as plt
import uproot
import pandas

model = "my_model_1"
input_root_file = 'root/spacepoints.root_output_invertedADC.root'
uproot_file = uproot.open(input_root_file)


tree = uproot_file['tsig_train']
featmap_file_name = "models/%s_featmap.txt"%(model)
featmap_file = open(featmap_file_name,'w')
featmap_list = []
featmap_counter = 0

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

    
bdt = xgb.Booster()
bdt.load_model("models/%s_model.json"%(model))


#t_test_sig = uproot_file['tsig_test'].arrays(featmap_list,library="pd",filter_branch=lambda b: b.dca.mid < 100000)
t_test_sig = uproot_file['tsig_test'].arrays(featmap_list,library="pd")
t_test_bkg = uproot_file['tbg_test'].arrays(featmap_list,library="pd")


xgb_test_sig  = xgb.DMatrix(t_test_sig,label=[1]*len(t_test_sig))
xgb_test_bkg  = xgb.DMatrix(t_test_bkg,label=[0]*len(t_test_bkg))



#########################################
#########################################
#### file output
#########################################
#########################################

output_root_file = "root/%s_BDT_scores_invertedADC.root"%(model)

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

