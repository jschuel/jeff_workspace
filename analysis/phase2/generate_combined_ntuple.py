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

day = sys.argv[1] #day of study
ring = sys.argv[2] #Choose LER, HER, or LUMI, as appropriate

#SKB_input = "~/data/phase2/SKB_data/07-%s-2018_interpolated.root"%(day) #input file with SKB parameters
SKB_input = "~/data/phase2/SKB_data/June_%s_SKB.root"%(day)
module_id = ["iiwi", "honu", "kohola", "nene", "tako", "humu", "palila", "elepaio"]
tpc_input = {} #MUST BE A DICTIONARY WITH MODULE_IDs AS KEYS
for module in module_id:
    tpc_input[module] = "~/data/phase2/background_studies/6-%s_%s/most_updated/%s_all.root"%(day,ring,module)
output_f = '~/data/phase2/combined_SKB_TPC_ntuples/June_%s_study.root'%(day) #output combined TPC + SKB analysis ntuple

make_ntuple(SKB_input, tpc_input, output_f) #generates the combined analysis ntuple from an SKB input file and a TPC input file. Parameters must be passed in this order
