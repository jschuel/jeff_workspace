from root_pandas import read_root
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ROOT
import array

#Build ROOT ntuple from combined dataframe
def make_ntuple(SKB_input, tpc_input, output_file):
    df = make_combined_dataframe(SKB_input,tpc_input)
    keys = [val for val in df.columns]
    output = ROOT.TFile(output_file, 'recreate')
    tout = ROOT.TTree('tout','tout')
    branches = {}
    data={}
    
    for key in keys:
        if df[key].dtype == "O": #Determines the size of an array in a dataframe to be pushed to ntuple
            npoints = len(df[key][0])
            data[key]=array.array('d',[0 for j in range(0,npoints)])
            branches[key]=tout.Branch("%s"%(key), data[key], "%s[%s]/D"%(key,npoints))
        else:
            data[key]=array.array('d',[0])
            branches[key]=tout.Branch("%s"%(key), data[key], "%s/D"%(key))

    for j in range(0,len(df)):
        for key in keys:
            if df[key].dtype == "O":
                npoints = len(df[key][0])
                for i in range(0,npoints):
                    data[key][i]=df[key][j][i]
            else:
                data[key][0]=df[key][j]
        tout.Fill()

    output.Write()
    output.Close()

#Merge TPC and SKB dataframes
def make_combined_dataframe(SKB_input,tpc_input):
    df_SKB = make_SKB_dataframe(SKB_input)
    SKB_ts = df_SKB['ts'].to_numpy()
    dfs = {} #dictionary of dataframes
    neutron_counts = {} #dictionary of neutron counts for each TPC to be merged with df_SKB
    total_counts = {}
    module_id = [key for key in tpc_input.keys()] #keys of tpc_input are defined to be module_ids
    for module in module_id:
        dfs[module] = make_TPC_dataframe(tpc_input, module)
        neutron_counts[module], total_counts[module] = merge_TPC_rates_with_SKB(dfs[module], SKB_ts)
        df_SKB['%s_neutrons'%(module)] = neutron_counts[module]
        df_SKB['%s_all_tracks'%(module)] = total_counts[module]
    return df_SKB
    
#Get TPC rates at SKB 1s time intervals
def merge_TPC_rates_with_SKB(df_TPC, ts_range):
    TPC_neutron_counts = []
    TPC_total_counts = []
    for i in range(0,len(ts_range)-1):
        if (ts_range[i+1]-ts_range[i]) <= 1.1:
            TPC_neutron_counts.append(len(df_TPC.loc[(df_TPC.timestamp_start > ts_range[i]) & (df_TPC.timestamp_start < ts_range[i+1]) & (df_TPC['neutron_flag'] == 1)].index))
            TPC_total_counts.append(len(df_TPC.loc[(df_TPC.timestamp_start > ts_range[i]) & (df_TPC.timestamp_start < ts_range[i+1])].index))
        else:
            TPC_neutron_counts.append(0)
            TPC_total_counts.append(0)
    TPC_neutron_counts.append(0)
    TPC_total_counts.append(0)
    return TPC_neutron_counts, TPC_total_counts

##Make TPC dataframe
def make_TPC_dataframe(tpc_input, module_id): #tpc_input is a dictionary of files with module_id's as keys
    df_TPC = read_root(tpc_input[module_id], "data")
    df_TPC['neutron_flag'] = [0 for i in range(0,len(df_TPC))]
    index = df_TPC.loc[(df_TPC.track_energy < (0.5*df_TPC.length-75)) & (df_TPC.track_energy > (0.04*df_TPC.length-65)) & (df_TPC.track_energy > 100) & (df_TPC.hitside_col_min == 0) & (df_TPC.hitside_col_max == 0) & (df_TPC.hitside_row_min == 0) & (df_TPC.hitside_row_max == 0)].index.to_numpy()
    df_TPC['neutron_flag'][index] = 1
    #df_TPC_neutron = df_TPC.iloc[df_TPC.loc[(df_TPC.track_energy < (0.5*df_TPC.length-75)) & (df_TPC.track_energy > (0.04*df_TPC.length-65)) & (df_TPC.track_energy > 100) & (df_TPC.hitside_col_min == 0) & (df_TPC.hitside_col_max == 0) & (df_TPC.hitside_row_min == 0) & (df_TPC.hitside_row_max == 0)].index] #dataframe for TPC nuclear recoils
    return df_TPC

##Make dataframe of SKB variables relevant for study
def make_SKB_dataframe(SKB_input):
    df_SKB = read_root(SKB_input)
    df_SKB = df_SKB.sort_values(by = ['ts']) #order by ascending timestamp
    df_SKB.index = [i for i in range(0,len(df_SKB))]
    return df_SKB
