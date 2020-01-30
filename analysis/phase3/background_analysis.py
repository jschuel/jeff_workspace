import root_pandas as rp
import pandas as pd
import ROOT
import math
import numpy as np
import matplotlib.pyplot as plt
from analysis_module import extract_variables
from analysis_module import make_heuristic_plots
from os import sys

month = sys.argv[1]
day = sys.argv[2]
ring = sys.argv[3]

input_file = '/Users/vahsengrouplaptop/data/phase3/combined_SKB_TPC_ntuples/%s_%s_%s_study.root'%(month, day, ring)
module_ids = ["iiwi", "humu", "nene", "tako", "palila", "elepaio"]
bin_width = 300

df_input = rp.read_root(input_file) #gives dataframe to pass into extract_variables

df = extract_variables(df_input, module_ids, ring, bin_width) #gives dataframe that's computed averages over bin_width slices of input df. Variables can be directly passed into combined heuristic for analysis

gr = make_heuristic_plots(df, module_ids)
