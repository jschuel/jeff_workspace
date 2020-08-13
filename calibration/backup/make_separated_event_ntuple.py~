from root_pandas import read_root
import pandas as pd
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sklearn.cluster as cluster
import array

from separate_events import create_ntuple #module uses all of the above
from os import sys

module_id = sys.argv[1]
run_num = sys.argv[2]
ring = sys.argv[3]

if ring=="HER":
    f_input = '~/tpc_tools_workspace/tpc_tools/output_ntuples/6-11_HER/initial_ntuples/%s_%s.root'%(module_id,run_num)
    f_output = '~/tpc_tools_workspace/tpc_tools/output_ntuples/6-11_HER/separated_ntuples/%s_%s.root'%(module_id, run_num)
else:
    f_input = '~/tpc_tools_workspace/tpc_tools/output_ntuples/6-12_LER/initial_ntuples/%s_%s.root'%(module_id,run_num)
    f_output = '~/tpc_tools_workspace/tpc_tools/output_ntuples/6-12_LER/separated_ntuples/%s_%s.root'%(module_id, run_num)

create_ntuple(f_input, f_output) #function from separate_events that creates the ntuple

