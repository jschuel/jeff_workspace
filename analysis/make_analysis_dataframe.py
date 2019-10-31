from root_pandas import read_root
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ROOT

#Merge TPC and SKB dataframes
def make_combined_dataframe(date, ring):
    if date == "11":
        study_indices = [i for i in range(73790,74920)] + [i for i in range(75138,76220)] + [i for i in range(76380,77520)] + [i for i in range(77650,78750)] +[i for i in range(78850,80670)] + [i for i in range(81570,84580)] + [i for i in range(84800,86390)] + [i for i in range(86550,87700)] + [i for i in range(87800,88650)] + [i for i in range(88900,89880)]
        Nb = [789 for i in range(73790, 74920)] + [789 for i in range(75138,76220)] + [789 for i in range(76380,77520)] + [789 for i in range(77650,78750)] + [789 for i in range(78850,80670)] + [1576 for i in range(81570,84580)] + [1576 for i in range(84800,86390)] + [1576 for i in range(86550,87700)] + [1576 for i in range(87800,88650)] + [1576 for i in range(88900,89880)]
        knob = [0 for i in range(73790, 74920)] + [1 for i in range(75138,76220)] + [-1 for i in range(76380,77520)] + [2 for i in range(77650,78750)] + [-2 for i in range(78850,80670)] + [0 for i in range(81570,84580)] + [0 for i in range(84800,86390)] + [1 for i in range(86550,87700)] + [2 for i in range(87800,88650)] + [0 for i in range(88900,89880)]

    if date == "12":
        study_indices = [i for i in range(126880,128170)] + [i for i in range(128490,129950)] + [i for i in range(130280,131600)] + [i for i in range(131910,132920)] + [i for i in range(133260,134120)] + [i for i in range(134370,135160)] + [i for i in range(138130,138920)] + [i for i in range(139100,142550)] + [i for i in range(142790,143600)] + [i for i in range(144080,145170)] + [i for i in range(145430,146020)] + [i for i in range(146160,146660)] + [i for i in range(146865,147245)] + [i for i in range(147540,149145)]
        Nb = [789 for i in range(126880,128170)] + [789 for i in range(128490,129950)] + [789 for i in range(130280,131600)] + [789 for i in range(131910,132920)] + [789 for i in range(133260,134120)] + [789 for i in range(134370,135160)] + [789 for i in range(138130,138920)] + [789 for i in range(139100,142550)] + [789 for i in range(142790,143600)] + [1576 for i in range(144080,145170)] + [1576 for i in range(145430,146020)] + [1576 for i in range(146160,146660)] + [1576 for i in range(146865,147245)] + [1576 for i in range(147540,149145)]
        knob = [0 for i in range(126880,128170)] + [0 for i in range(128490,129950)] + [0 for i in range(130280,131600)] + [0 for i in range(131910,132920)] + [1 for i in range(133260,134120)] + [-1 for i in range(134370,135160)] + [1.9 for i in range(138130,138920)] + [1.9 for i in range(139100,142550)] + [1.9 for i in range(142790,143600)] + [0 for i in range(144080,145170)] + [1 for i in range(145430,146020)] + [1.9 for i in range(146160,146660)] + [-1 for i in range(146865,147245)] + [-1.9 for i in range(147540,149145)]
    df_SKB = make_SKB_dataframe(date)
    df_study = df_SKB.iloc[study_indices]
    df_study['Nb'] = Nb
    df_study['knob'] = knob
    SKB_ts = df_study['ts'].to_numpy()
    df_iiwi = make_TPC_dataframe(date, ring, "iiwi")
    df_honu = make_TPC_dataframe(date, ring, "honu")
    df_kohola = make_TPC_dataframe(date, ring, "kohola")
    df_nene = make_TPC_dataframe(date, ring, "nene")
    df_tako = make_TPC_dataframe(date, ring, "tako")
    df_humu = make_TPC_dataframe(date, ring, "humu")
    df_palila = make_TPC_dataframe(date, ring, "palila")
    df_elepaio = make_TPC_dataframe(date, ring, "elepaio")
    iiwi_neutron_counts = merge_TPC_rates_with_SKB(df_iiwi, SKB_ts)
    df_study['iiwi_neutrons'] = iiwi_neutron_counts
    honu_neutron_counts = merge_TPC_rates_with_SKB(df_honu, SKB_ts)
    df_study['honu_neutrons'] = honu_neutron_counts
    kohola_neutron_counts = merge_TPC_rates_with_SKB(df_kohola, SKB_ts)
    df_study['kohola_neutrons'] = kohola_neutron_counts
    nene_neutron_counts = merge_TPC_rates_with_SKB(df_nene, SKB_ts)
    df_study['nene_neutrons'] = nene_neutron_counts
    tako_neutron_counts = merge_TPC_rates_with_SKB(df_tako, SKB_ts)
    df_study['tako_neutrons'] = tako_neutron_counts
    humu_neutron_counts = merge_TPC_rates_with_SKB(df_humu, SKB_ts)
    df_study['humu_neutrons'] = humu_neutron_counts
    palila_neutron_counts = merge_TPC_rates_with_SKB(df_palila, SKB_ts)
    df_study['palila_neutrons'] = palila_neutron_counts
    elepaio_neutron_counts = merge_TPC_rates_with_SKB(df_elepaio, SKB_ts)
    df_study['elepaio_neutrons'] = elepaio_neutron_counts
    return df_study
    
#Get TPC rates at SKB 1s time intervals
def merge_TPC_rates_with_SKB(df_TPC, ts_range):    
    TPC_neutron_counts = []
    for i in range(0,len(ts_range)-1):
        if (ts_range[i+1]-ts_range[i]) <= 2:
            TPC_neutron_counts.append(len(df_TPC.loc[(df_TPC.ts > ts_range[i]) & (df_TPC.ts < ts_range[i+1])].index))
        else:
            TPC_neutron_counts.append(0)
    TPC_neutron_counts.append(0)
    return TPC_neutron_counts

##Make TPC dataframe
def make_TPC_dataframe(date, ring, module_id):
    df_TPC = read_root("/Users/vahsengrouplaptop/data/phase2/background_studies/6-%s_%s/tpc_tools_processed/separated_ntuples/whole_study_separated_%s.root"%(date, ring, module_id),"tracks")
    df_TPC_neutron = df_TPC.iloc[df_TPC.loc[(df_TPC.recoil_energy < (0.25*df_TPC.length-75)) & (df_TPC.recoil_energy > (0.0411764*df_TPC.length-64.688)) & (df_TPC.recoil_energy > 100)].index] #dataframe for TPC nuclear recoils
    return df_TPC_neutron

##Make dataframe of SKB variables relevant for study
def make_SKB_dataframe(date):
    df_SKB = read_root("/Users/vahsengrouplaptop/data/phase2/SKB_data/SKB_processed_2018-06-%s_all_pressures.root"%(date))
    df_SKB = df_SKB.sort_values(by = ['ts']) #order by ascending timestamp
    df_SKB.index = [i for i in range(0,len(df_SKB))]
    return df_SKB

