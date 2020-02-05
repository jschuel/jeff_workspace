from root_pandas import read_root
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import ROOT
import array

#Build ROOT ntuple from combined dataframe
def make_ntuple(SKB_input, tpc_input, study_indices, output_file):
    df = make_combined_dataframe(SKB_input,tpc_input,study_indices)
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
def make_combined_dataframe(SKB_input,tpc_input,study_indices):
    df_SKB = make_SKB_dataframe(SKB_input)
    SKB_ts = df_SKB['ts'].to_numpy()
    dfs = {} #dictionary of dataframes
    neutron_counts = {} #dictionary of neutron counts for each TPC to be merged with df_SKB
    module_id = [key for key in tpc_input.keys()] #keys of tpc_input are defined to be module_ids
    for module in module_id:
        dfs[module] = make_TPC_dataframe(tpc_input, module)
        neutron_counts[module] = merge_TPC_rates_with_SKB(dfs[module], SKB_ts)
        df_SKB['%s_neutrons'%(module)] = neutron_counts[module]
    df_SKB['Storage_Flag'] = [0 for i in range(0,len(df_SKB))]
    df_SKB['Storage_Flag'][study_indices] = 1
    return df_SKB
    
#Get TPC rates at SKB 1s time intervals
def merge_TPC_rates_with_SKB(df_TPC, ts_range):
    #if month == "May":
    #    if side == "BWD":
    #        df_TPC.ts = df_TPC.ts + 213 #to account for BWD TPCs being 213 seconds behind NTP. Fixxed for Autumn 2019 runs and beyond
    TPC_neutron_counts = []
    for i in range(0,len(ts_range)-1):
        if (ts_range[i+1]-ts_range[i]) <= 1.1:
            TPC_neutron_counts.append(len(df_TPC.loc[(df_TPC.ts > ts_range[i]) & (df_TPC.ts < ts_range[i+1])].index))
        else:
            TPC_neutron_counts.append(0)
    TPC_neutron_counts.append(0)
    return TPC_neutron_counts

##Make TPC dataframe
def make_TPC_dataframe(tpc_input, module_id): #tpc_input is a dictionary of files with module_id's as keys
    df_TPC = read_root(tpc_input[module_id], "tracks")
    df_TPC['recoil_energy'] = df_TPC['recoil_energy']/1000
    df_TPC_neutron = df_TPC.iloc[df_TPC.loc[(df_TPC.recoil_energy < (0.5*df_TPC.length-75)) & (df_TPC.recoil_energy > (0.0411764*df_TPC.length-64.688)) & (df_TPC.recoil_energy > 100) & (df_TPC.hitside_top == 0) & (df_TPC.hitside_bottom == 0) & (df_TPC.hitside_source == 0) & (df_TPC.hitside_antisource == 0)].index] #dataframe for TPC nuclear recoils
    return df_TPC_neutron

##Make dataframe of SKB variables relevant for study
def make_SKB_dataframe(SKB_input):
    df_SKB = read_root(SKB_input)
    df_SKB = df_SKB.drop(columns=['HE3', 'TPC'])
    df_SKB = df_SKB.sort_values(by = ['ts']) #order by ascending timestamp
    df_SKB.index = [i for i in range(0,len(df_SKB))]
    return df_SKB

#Study indices is a parameter passed into the ntuple builder. It sets the storage flag for the combined ntuple.
#User can define boolean expressions for appropriate background studies. Should come up with a more sophisticated way to do this in the future
def get_study_indices(month,day,ring):
    if month == "May": #Study indices for generating storage flag
        if day == "11":
            study_indices = [i for i in range(10586,12258)] + [i for i in range(12972,13828)] + [i for i in range(14443,14764)] + [i for i in range(15180,15508)] + [i for i in range(15658,16075)] + [i for i in range(16731,17178)] + [i for i in range(17667,18142)] + [i for i in range(18989,19496)] + [i for i in range(19923,20509)] + [i for i in range(21043,22085)] + [i for i in range(22465,25085)] + [i for i in range(25815,27231)] + [i for i in range(27531,27988)] + [i for i in range(29311,30599)]
        
        if day == "12":
            study_indices = [i for i in range(1598,2636)] + [i for i in range(2756,4425)] + [i for i in range(4655,5884)] + [i for i in range(7269,8552)] + [i for i in range(8672,10182)] + [i for i in range(10502,11544)] + [i for i in range(12381,12985)] + [i for i in range(13121,13911)] + [i for i in range(14086,14675)] + [i for i in range(15711,16753)] + [i for i in range(17013,17997)] + [i for i in range(18076,18951)] + [i for i in range(19896,20731)] + [i for i in range(20811,21645)] + [i for i in range(21805,22739)]

        if day == "14": #omitted indices (10838-11981) these have changing YaECK but they're in too short of bursts for TPC rates
            study_indices = [i for i in range(1202,2736)] + [i for i in range(2976,3729)] + [i for i in range(4317,5757)] + [i for i in range(6241,6514)] + [i for i in range(6704,7726)] + [i for i in range(8488,8651)] + [i for i in range(8831,10092)] + [i for i in range(12290,12974)] + [i for i in range(13154,13472)] + [i for i in range(13552,14069)] + [i for i in range(14179,14764)] + [i for i in range(14914,15410)]

    if month == "Dec":
        if ring == "LER":
            study_indices = [i for i in range(280,2190)] + [i for i in range(3920,5425)] + [i for i in range(6850,7665)] + [i for i in range(8345,8675)] + [i for i in range(9090,9420)] + [i for i in range(9820,10180)] + [i for i in range(10600,11040)] + [i for i in range(11600,13180)] + [i for i in range(13440,13810)] + [i for i in range(13960,14440)] + [i for i in range(14975,15185)] + [i for i in range(15360,15580)] + [i for i in range(15745,16050)] + [i for i in range(16775,18110)] + [i for i in range(18540,19610)] + [i for i in range(19980,21150)]
        if ring == "HER":
            study_indices = [i for i in range(580,2510)] + [i for i in range(3465,4930)] + [i for i in range(5560,7390)] + [i for i in range(8000,10030)] + [i for i in range(10375,12115)] + [i for i in range(14100,16075)] + [i for i in range(16700,18670)] + [i for i in range(22680,25354)]
        if ring =="LUMI":
            study_indices = [i for i in range(0,23405)]
    return study_indices
