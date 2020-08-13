'''JTS: 01-29-2020: Uses the ntuple_maker_module to pass in a scan .h5 file and generate an preliminary ntuple that will be passed to tpc_tools

'''

import ROOT
import array
import numpy as np
import pandas as pd
from ntuple_maker_module import create_ntuple
import sys


month = sys.argv[1]
day = sys.argv[2]
speed = sys.argv[3]

#sides = ["FWD", "BWD"]
sides = ["BWD"]
for side in sides:
    if side == "FWD":
        TPCs = ["nene","humu","iiwi"]
        #run_num = [i for i in range(24740,24749)] + [i for i in range(24827,24890)] #for dec 7th study
        run_num = [i for i in range(34745,34820)]
    if side == "BWD":
        #TPCs = ["elepaio", "palila", "tako"]
        TPCs = ["elepaio"]
        run_num = [346]
        #run_num = [343,344,345,346] #for dec 7th study
        #run_num = [348,349,350] #for dec 8th study
    for tpc in TPCs:
        for run in run_num:
            f_input = '~/data/phase3/phase3_background_h5/%s_%s/%s_%s_stop_mode_ext_trigger_scan_interpreted.h5'%(month, day, run, tpc) 
            #f_output = '~/data/phase3/phase3_background_root/%s_%s/%s_%s.root'%(month, day, run, tpc)
            f_output = 'output/%s/%s_%s.root'%(speed,tpc,run)
            create_ntuple(f_input, f_output, speed)
