#!/usr/bin/python3

import sys
import os

def preselection(f,sample_dir,tag):
    print("Python wrapper:")
    print('root -l -b -q macro/preselection.C\'("%s","%s","%s")\''%(f,sample_dir,tag))
    output_dir=sample_dir.replace('raw_sample','preselection')
    os.system('mkdir -p %s/%s'%(output_dir,tag))
    os.system('mkdir -p %s/%s'%(output_dir.replace('samples/','img/'),tag))
    os.system('root -l -b -q macro/preselection.C\'("%s","%s","%s")\''%(sample_dir,f,tag))

'''
takes 2 arguments, sample and tag then loops over all files in the sample
if you have a sample MillichargeGamma3D_2406 run as

$ python3 run_preselection.py samples/MillichargeGamma3D_2406 mypreselection_240626

if you want to run detvars for this sample, just use that directory

$ python3 run_preselection.py samples/MillichargeGamma3D_2406/detvars mypreselection_240626
'''

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

if len(sys.argv) > 2:
    sample_dir = sys.argv[1]
    tag        = sys.argv[2]
    print("running sample %s for all files:"%(sample_dir))
    counter = 1
    for f in os.listdir("%s"%(sample_dir)):
        if not f.endswith(".root"):
            print("ERROR: sample dir should be where the raw .root files are located")
            break
        counter = counter+1
        preselection(f,sample_dir,tag)
        if counter > 1:
            break
elif len(sys.argv) == 2:
    sample_dir = sys.argv[1]
    if sample_dir.endswith('.root'):
        print("running single file %s"%(sample_dir))
        preselection(sample_dir)
    else:
        print("Need two arguments: sample, tag. One argument only for single file ending in .root")
else:
    print("sample name as first argument, preselection tag as second argument")
