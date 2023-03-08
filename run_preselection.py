#!/usr/bin/python3

import sys
import os

def preselection(f):
    os.system('root -l -b -q macro/preselection.C\'("%s")\''%f)
    

# three run modes, single file, masses, systematics

masses = [ 100, 150, 200, 300, 350, 400 ]
detvars = [ 'CV',
            'wiremod_x',
            'wiremod_yz',
            'wiremod_anglexz',
            'wiremod_angleyz',
            'wiremod_dEdx',
            'ly_down',
            'ly_rayleigh',
            'ly_atten'
           ]


if len(sys.argv) > 1:
    arg = sys.argv[1]
    print(arg)

    if arg.endswith('.root'):
        preselection(arg)
        exit()
    elif "detvar" in arg:
        for dv in detvars:
            preselection("systematics/%s/spacepoints_%s.root"%(arg,dv))
    elif "mev" in arg:
        for m in masses:
            preselection("root/%s/spacepoints_%smev.root"%(arg,m))
    else:
        preselection("root/%s/spacepoints_%s.root"%(arg,arg))
