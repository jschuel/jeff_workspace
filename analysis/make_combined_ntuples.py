###Script creates ROOT ntuples that combine SKB and TPC data for background studies of interest. make_dataframe module makes a pandas dataframe with relevant SKB variables as well as corresponding TPC rates for all TPCs. This script then uses root_pandas to save the dataframe as an ntuple

from root_pandas import to_root
from make_analysis_dataframe import *

dates = ["11","12"] #11 is HER, 12 is LER study

for i in dates:
    if i == "11":
        df = make_combined_dataframe(i, "HER")
        df.to_root('/Users/vahsengrouplaptop/data/phase2/combined_SKB_TPC_ntuples/June_11_HER.root', key='data')
    if i == "12":
        df = make_combined_dataframe(i, "LER")
        df.to_root('/Users/vahsengrouplaptop/data/phase2/combined_SKB_TPC_ntuples/June_12_LER.root', key='data')

