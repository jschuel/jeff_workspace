import root_pandas as rp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_data(input_file):
    df = rp.read_root(input_file, key = 'truthB')
    return df

def add_delta_t(df): #input dataframe
    dt = np.array([np.abs(df.loc[df['__event__']==event]['mcFlightTime'].diff().iloc[1]) for event in df['__event__']])
    df['dt'] = dt*1000
    return df

def determine_same_flavor(df):
    pdg_sum = np.array([np.abs(df.loc[df['__event__']==event]['mcPDG'].sum()) for event in df['__event__']])
    df['pdg_sum'] = pdg_sum
    df['OF'] = 0
    df['SF'] = 0
    index = df.loc[df['pdg_sum'] == 0].index.to_numpy()
    df['OF'][index] = 1
    df['SF'][df.loc[df['OF'] == 0].index.to_numpy()] = 1
    df = df.drop(columns = ['pdg_sum'])
    return df

def reduce_and_bin_data(df):
    df = add_delta_t(df)
    #df = df.rename(columns={'MCDeltaT':'dt'})
    df = determine_same_flavor(df)
    df['bin'] = 0
    df_red = df[['OF', 'SF', 'dt', 'bin']]
    df_red['dt'] = np.abs(df_red['dt'])
    df_red = df_red.loc[df_red['dt'].duplicated()==False]
    df_red['bin'][df_red.loc[(df_red['dt']>=0.5) & (df_red['dt']< 1)].index.to_numpy()]=1 #Bin according to Go's convention
    df_red['bin'][df_red.loc[(df_red['dt']>=1) & (df_red['dt']< 2)].index.to_numpy()]=2
    df_red['bin'][df_red.loc[(df_red['dt']>=2) & (df_red['dt']< 3)].index.to_numpy()]=3
    df_red['bin'][df_red.loc[(df_red['dt']>=3) & (df_red['dt']< 4)].index.to_numpy()]=4
    df_red['bin'][df_red.loc[(df_red['dt']>=4) & (df_red['dt']< 5)].index.to_numpy()]=5
    df_red['bin'][df_red.loc[(df_red['dt']>=5) & (df_red['dt']< 6)].index.to_numpy()]=6
    df_red['bin'][df_red.loc[(df_red['dt']>=6) & (df_red['dt']< 7)].index.to_numpy()]=7
    df_red['bin'][df_red.loc[(df_red['dt']>=7) & (df_red['dt']< 9)].index.to_numpy()]=8
    df_red['bin'][df_red.loc[(df_red['dt']>=9) & (df_red['dt']< 13)].index.to_numpy()]=9
    df_red['bin'][df_red.loc[(df_red['dt']>=13) & (df_red['dt']< 20)].index.to_numpy()]=10
    df_red['asymmetry'] = 0
    asymmetry = []
    for i in range(0,11):
        asymmetry.append((df_red.loc[df_red['bin']==i]['OF'].sum()-df_red.loc[df_red['bin']==i]['SF'].sum())/(df_red.loc[df_red['bin']==i]['OF'].sum()+df_red.loc[df_red['bin']==i]['SF'].sum()))
        index = df_red.loc[df_red['bin'] == i].index.to_numpy()
        df_red['asymmetry'][index] = asymmetry[i]
    print(asymmetry)
    return df_red

data = get_data("~/truth_signal.root")
df_red = reduce_and_bin_data(data)
df_red.to_root("digested_truth.root", key = 'B')

    
