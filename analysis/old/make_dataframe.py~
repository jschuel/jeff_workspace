from root_pandas import read_root
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ROOT

#Merge TPC and SKB dataframes
def make_combined_dataframe(date, ring):
    if date == "11":
        study_indices = [i for i in range(10586,12258)] + [i for i in range(12972,13828)] + [i for i in range(14443,14764)] + [i for i in range(15180,15508)] + [i for i in range(15658,16075)] + [i for i in range(16731,17178)] + [i for i in range(17667,18142)] + [i for i in range(18989,19496)] + [i for i in range(19923,20509)] + [i for i in range(21043,22085)] + [i for i in range(22465,25085)] + [i for i in range(25815,27231)] + [i for i in range(27531,27988)] + [i for i in range(29311,30599)]

    if date == "12":
        study_indices = [i for i in range(1598,2636)] + [i for i in range(2756,4425)] + [i for i in range(4655,5884)] + [i for i in range(7269,8552)] + [i for i in range(8672,10182)] + [i for i in range(10502,11544)] + [i for i in range(12381,12985)] + [i for i in range(13121,13911)] + [i for i in range(14086,14675)] + [i for i in range(15711,16753)] + [i for i in range(17013,17997)] + [i for i in range(18076,18951)] + [i for i in range(19896,20731)] + [i for i in range(20811,21645)] + [i for i in range(21805,22739)]

    if date == "14": #omitted indices (10838-11981) these have changing YaECK but they're in too short of bursts for TPC rates
        study_indices = [i for i in range(1202,2736)] + [i for i in range(2976,3729)] + [i for i in range(4317,5757)] + [i for i in range(6241,6514)] + [i for i in range(6704,7726)] + [i for i in range(8488,8651)] + [i for i in range(8831,10092)] + [i for i in range(12290,12974)] + [i for i in range(13154,13472)] + [i for i in range(13552,14069)] + [i for i in range(14179,14764)] + [i for i in range(14914,15410)]

    df_SKB = make_SKB_dataframe(date)
    df_study = df_SKB.iloc[study_indices]
    SKB_ts = df_study['ts'].to_numpy()
    df_humu = make_TPC_dataframe(date, ring, "humu")
    df_nene = make_TPC_dataframe(date, ring, "nene")
    df_tako = make_TPC_dataframe(date, ring, "tako")
    df_palila = make_TPC_dataframe(date, ring, "palila")
    df_elepaio = make_TPC_dataframe(date, ring, "elepaio")
    humu_neutron_counts = merge_TPC_rates_with_SKB(df_humu, SKB_ts, "FWD")
    df_study['humu_neutrons'] = humu_neutron_counts
    nene_neutron_counts = merge_TPC_rates_with_SKB(df_nene, SKB_ts, "FWD")
    df_study['nene_neutrons'] = nene_neutron_counts
    tako_neutron_counts = merge_TPC_rates_with_SKB(df_tako, SKB_ts, "BWD")
    df_study['tako_neutrons'] = tako_neutron_counts
    palila_neutron_counts = merge_TPC_rates_with_SKB(df_palila, SKB_ts, "BWD")
    df_study['palila_neutrons'] = palila_neutron_counts
    elepaio_neutron_counts = merge_TPC_rates_with_SKB(df_elepaio, SKB_ts, "BWD")
    df_study['elepaio_neutrons'] = elepaio_neutron_counts
    return df_study
    
#Get TPC rates at SKB 1s time intervals
def merge_TPC_rates_with_SKB(df_TPC, ts_range, side):
    if side == "BWD":
        df_TPC.ts = df_TPC.ts + 213 #to account for BWD TPCs being 213 seconds behind NTP    
    TPC_neutron_counts = []
    for i in range(0,len(ts_range)-1):
        if (ts_range[i+1]-ts_range[i]) <= 1.1:
            TPC_neutron_counts.append(len(df_TPC.loc[(df_TPC.ts > ts_range[i]) & (df_TPC.ts < ts_range[i+1])].index))
        else:
            TPC_neutron_counts.append(0)
    TPC_neutron_counts.append(0)
    return TPC_neutron_counts

##Make TPC dataframe
def make_TPC_dataframe(date, ring, module_id):
    df_TPC = read_root("/Users/vahsengrouplaptop/data/phase3/phase3_background_root/tpc_tools/May_%s_%s_%s_phase3.root"%(date, ring, module_id),"tracks")
    df_TPC_neutron = df_TPC.iloc[df_TPC.loc[(df_TPC.recoil_energy < (0.25*df_TPC.length-75)) & (df_TPC.recoil_energy > (0.0411764*df_TPC.length-64.688)) & (df_TPC.recoil_energy > 100)].index] #dataframe for TPC nuclear recoils
    return df_TPC_neutron

##Make dataframe of SKB variables relevant for study
def make_SKB_dataframe(date):
    if date == "11" or date == "14":
        df_SKB = read_root("/Users/vahsengrouplaptop/data/phase3/PVM/May_%s_PVM.root"%(date), columns=['ts', 'SKB_BMLDCCT_CURRENT', 'SKB_BMLXRM_BEAM_SIGMAY', 'SKB_VALCCG_LER_PRES_AVG', 'SKB_CGLINJ_BKSEL_NOB_SET'])
    else:
        df_SKB = read_root("/Users/vahsengrouplaptop/data/phase3/PVM/May_%s_PVM.root"%(date), columns=['ts', 'SKB_BMHDCCT_CURRENT', 'SKB_BMHXRM_BEAM_SIGMAY', 'SKB_VAHCCG_HER_PRES_AVG', 'SKB_CGHINJ_BKSEL_NOB_SET'])
    df_SKB.columns = ['ts', 'I', 'sigy', 'P', 'Nb']
    df_SKB = df_SKB.sort_values(by = ['ts']) #order by ascending timestamp
    df_SKB.index = [i for i in range(0,len(df_SKB))]
    return df_SKB

