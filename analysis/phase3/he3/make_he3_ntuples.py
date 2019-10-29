###Script creates ROOT ntuples that combine SKB and TPC data for background studies of interest. make_dataframe module makes a pandas dataframe with relevant SKB variables as well as corresponding TPC rates for all TPCs. This script then uses root_pandas to save the dataframe as an ntuple

from root_pandas import to_root
from merge_he3 import *

dates = ["11","12","14"] #11 and 14 are LER, 12 is HER study

for i in dates:
    if i == "11":
        df = make_combined_dataframe(i, "LER")
        df.to_root('/Users/vahsengrouplaptop/data/phase3/combined_SKB_TPC_ntuples/May_11_LER_He3.root', key='data')
    if i == "12":
        df = make_combined_dataframe(i, "HER")
        df.to_root('/Users/vahsengrouplaptop/data/phase3/combined_SKB_TPC_ntuples/May_12_HER_He3.root', key='data')
    if i == "14":
        df = make_combined_dataframe(i, "LER")
        df.to_root('/Users/vahsengrouplaptop/data/phase3/combined_SKB_TPC_ntuples/May_14_LER_He3.root', key='data')
