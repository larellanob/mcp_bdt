import uproot
import pandas
from ROOT import TH1F
from ROOT import TFile

location = 'outputs'
sig_files = [f'{location}/checkout_train_5MeV.root',] \
    + [f'{location}/checkout_train{i}_5MeV.root' for i in range(2,21)]
bkg_files = [f'{location}/checkout_prodgenie_numi_nu_overlay_run1.root',
             f'{location}/checkout_data_extnumi_run1.root']

# concatenate files
# leave space to filter by "used in training"
# need to grab pot from input files
# run make_bdt_histograms.cxx

sig_bdt = pandas.concat([ uproot.open(s)['tree']['bdt_result'].arrays(library='pd') for s in sig_files],ignore_index=True)
bkg_bdt = [ uproot.open(b)['tree']['bdt_result'].arrays(library='pd') for b in bkg_files ]
#for s in sig_files:
#    f = uproot.open(s)
#    bdt=f['tree']['bdt_result'].arrays(library='pd')


# grab pot and bkg nevents information from the input files
sig_pot = 0.0
bkg_pot = 0.0
bkg_eve = 0
in_location = 'inputs'
in_sig_files = [f'{in_location}/checkout_train_5MeV.root',] \
    + [f'{in_location}/checkout_train{i}_5MeV.root' for i in range(2,21)]
in_bkg_files = [f'{in_location}/checkout_prodgenie_numi_nu_overlay_run1.root',
             f'{in_location}/checkout_data_extnumi_run1.root']
for s in in_sig_files:
    f = uproot.open(s)
    pot_array = f['wcpselection/T_pot']['pot_tor875'].arrays(library='pd')
    #print(pot_array)
    sig_pot += pot_array.sum()
f = uproot.open(in_bkg_files[0])
bkg_pot_array = f['wcpselection/T_pot']['pot_tor875'].arrays(library='pd')
bkg_pot = bkg_pot_array.sum()


f = uproot.open(in_bkg_files[1])
bkg_eve = f['wcpselection/T_eval'].num_entries


print(sig_pot)
print(bkg_pot)
print(bkg_eve)

# i am going to export multiple backgrounds, these will have to be
# added as a list into the sensitivity thing, maybe make a class that
# has a member a list of backgrounds. i don't think you can make a
# proper list of histograms in root, or maybe you can, as an
# std::vector<TH1F*> or something, i'm going to try that but if it
# doesn't work i'll just have to use the names: bkg0 (overlay) bkg1
# (neutrino) bkg2 (whatever), etc. and the limits macro will have to
# recognize that somehow

# save to root file
# 1. sig histogram of bdt distribution
#    you know what, let's make it generalized once and for all, make a
#    vector or array of signals as well
# 2. array or std::vector of bkg histograms of bdt distribution
# 3. pot ttree with two branches: signal pots, bkg pots, bkg number of
#    events, and flag specifying which of these the bkg should use

# list of tuples (array,name to be saved)
sig_list = [ (sig_bdt,'test_sig0') ]
print("len of bkg_bdt",len(bkg_bdt))
bkg_list = [ (bkg_bdt[0],'test_bkg0'), (bkg_bdt[1],'test_bkg1')]

samples_list = sig_list+bkg_list
pot_list = [ ([sig_pot],'pot_sig0'), ([bkg_pot],'pot_bkg0'), ([bkg_eve],'evt_bkg1') ]




out_filename = 'test.root'

with uproot.recreate(out_filename) as f2:
    # first (this) loop is for trees
    branches_pot = {}
    for i in range(0,len(samples_list)):
        # if you need multiple branches, use another loop
        branches = {}
        print("type of samples:",type(samples_list[i][0]))
        branches['bdt'] = samples_list[i][0]
        f2[samples_list[i][1]] = branches
        print("type of pot:",type(pot_list[i][0]))
        branches_pot[pot_list[i][1]] = pot_list[i][0]
    f2["total_pot"] = branches_pot
    #f2["test_histo"] = h

rootfile = TFile(out_filename, "update")
h = TH1F("myHist", "myTitle", 64, -4, 4)
h.FillRandom("gaus",10000)
h1 = TH1F("sig0_hist", "Test BDT sig;Score;Entries", 40, 0, 1)
h2 = TH1F("bkg0_hist", "Test BDT bkg0;Score;Entries", 40, 0, 1)
h3 = TH1F("bkg1_hist", "Test BDT bkg1;Score;Entries", 40, 0, 1)
hist_list = [h1,h2,h3]
counter = 0
for h,s in zip(hist_list,samples_list):
    for i in s[0]['bdt_result']:
        h.Fill(i)

rootfile.Write()
