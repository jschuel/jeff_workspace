from root_pandas import read_root
import pandas as pd
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sklearn.cluster as cluster
import array

from separate_events import create_ntuple #module uses all of the above                                                                           
from os import sys

side = sys.argv[1]
if side == "FWD":
    TPCs = ["nene","humu"]
    run_num=[2940,2941,2942,2943,2944,2945,2946,2947,2963,2964,2965,2966,2967,2968,2969,2970,3012,3013,3014,3015,3016,3017,3018]
if side == "BWD":
    TPCs = ["elepaio", "palila", "tako"]
    run_num=[1554,1555,1556,1557,1558,1559,1560,1578,1579,1580,1581,1582,1583,1627,1628,1629,1630,1631,1632]

for tpc in TPCs:
    for run in run_num:
        f_input = '/Users/vahsengrouplaptop/data/phase3/phase3_background_root/%s_%s.root'%(run, tpc)
        f_output = '/Users/vahsengrouplaptop/data/phase3/phase3_background_root/%s_%s_separated_events.root'%(run, tpc)
        create_ntuple(f_input, f_output)



