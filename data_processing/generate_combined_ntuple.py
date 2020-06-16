'''JTS 01-29-2020
This script generates a combined TPC + SKB ntuple. Ideally this will be used to combine TPC data with pre-generated global ntuples.
'''

from root_pandas import read_root
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ROOT
import array
from os import sys
from make_combined_ntuple_module import make_ntuple #custom module for creating these ntuples
from make_combined_ntuple_module import get_study_indices

SKB_input = "/Users/vahsengrouplaptop/data/phase3/spring_2020/05-09-20/SKB+LUMI_archiver_with_flag.root"
module_id = ["tako", "palila", "elepaio", "iiwi", "nene", "humu"]
tpc_input = {}
for module in module_id:
    tpc_input[module] = "/Users/vahsengrouplaptop/data/phase3/spring_2020/05-09-20/tpc_root_files/%s_all_new.root"%(module)
study_indices = [i for i in range(0,34000)]
output_f = '/Users/vahsengrouplaptop/data/phase3/spring_2020/05-09-20/combined_ntuples/05-09_whole_study_newest.root'
make_ntuple(SKB_input, tpc_input, study_indices, output_f)
'''
month = sys.argv[1] #month of study
day = sys.argv[2] #day of study
ring = sys.argv[3] #Choose LER, HER, or LUMI, as appropriate

SKB_input = "/Users/vahsengrouplaptop/data/phase3/PVM/Dec_%s_%s_updated.root"%(day, ring) #input file with SKB parameters
if ring == "LUMI":
    module_id = ["tako", "palila", "elepaio", "iiwi", "nene", "humu"]
else:
    module_id = ["tako", "palila", "elepaio", "iiwi", "nene", "humu"] #Edit these for which TPCs you'd like included in the ntuple
tpc_input = {} #MUST BE A DICTIONARY WITH MODULE_IDs AS KEYS
for module in module_id:
    tpc_input[module] = "/Users/vahsengrouplaptop/data/phase3/phase3_background_root/tpc_tools/%s_%s_%s_phase3.root"%(month, day, module)
study_indices = get_study_indices(month, day, ring) #set indices for storage fills in a study. Add to this function in make_combined_ntuple.py for each new study
output_f = '/Users/vahsengrouplaptop/data/phase3/combined_SKB_TPC_ntuples/Dec_%s_%s_study_new.root'%(day, ring) #output combined TPC + SKB analysis ntuple

make_ntuple(SKB_input, tpc_input, study_indices, output_f) #generates the combined analysis ntuple from an SKB input file and a TPC input file. Parameters must be passed in this order
'''
