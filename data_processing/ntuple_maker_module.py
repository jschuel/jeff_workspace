''' 01/29/2020 JTS: Takes an input h5 file and converts it into a root file. 
After the root file is made, it can be passed into TPC_tools for final processing.
Branches generated with this function are:

event_number
timestamp_start
timestamp_stop
npoints
column
row
BCID
tot
sum_tot
num_clusters

The num_clusters branch gives the number of clusters in the event before it
was separated into individual clusters. This step is done because some events
contain multiple tracks and should be separated into individual tracks for PID.
'''
import pandas as pd
import numpy as np
import array
import matplotlib.pyplot as plt
import sklearn.cluster as cluster
import ROOT

def create_ntuple(f_input, f_output, speed):
    if speed == "fast":
        input_df = quick_h5_extract(f_input)
        maximum = 10000
    elif speed == "slow":
        input_df = extract_parameters_from_h5_file(f_input)
        maximum = 100000
    else:
        print("Please enter fast or slow")
    df = separate_events(input_df) #dataframe with multi-track events separated into single events. Call df for short
    #order: event_number, column, row, tot, BCID, timestamp, npoints, sum_tot, num_clusters
    output = ROOT.TFile(f_output, 'recreate')
    tree = ROOT.TTree('data', 'data')
    root_evt = array.array('i',[0])
    root_ts = array.array('d',[0])
    root_npoints = array.array('i',[0])
    root_sum_tot = array.array('i',[0])
    root_num_clusters = array.array('i',[0])
    root_col = array.array('i' , [0 for i in range(0,maximum)])
    root_row = array.array('i' , [0 for i in range(0,maximum)])
    root_tot = array.array('i' , [0 for i in range(0,maximum)])
    root_BCID = array.array('i' , [0 for i in range(0,maximum)])
    tree.Branch("event_number", root_evt, "event_number/I")
    tree.Branch("timestamp", root_ts , "timestamp/D")
    tree.Branch("npoints", root_npoints , "npoints/I")
    tree.Branch("column", root_col, "column[npoints]/I")
    tree.Branch("row", root_row, "row[npoints]/I")
    tree.Branch("BCID", root_BCID, "BCID[npoints]/I")
    tree.Branch("tot", root_tot, "tot[npoints]/I")
    tree.Branch("sum_tot", root_sum_tot, "sum_tot/I")
    tree.Branch("num_clusters", root_num_clusters, "num_clusters/I")

    for j in range(0,len(df)):
        root_evt[0] = df['event_number'].to_numpy()[j].astype(int)
        root_ts[0] = df['timestamp'].to_numpy()[j]
        root_npoints[0] = df['npoints'].to_numpy()[j].astype(int)
        root_sum_tot[0] = df['sum_tot'].to_numpy()[j].astype(int)
        root_num_clusters[0] = df['num_clusters'].to_numpy()[j].astype(int)
        for i in range(0,root_npoints[0]):
            root_col[i] = df['column'][j][i]-1
            root_row[i] = df['row'][j][i]-1
            root_tot[i] = df['tot'][j][i]
            root_BCID[i] = df['BCID'][j][i]
        tree.Fill()

    tree.Print()
    output.Write()
    output.Close()

def separate_events(dataframe): #takes read in dataframe and separates multitrack events using DBSCAN with epsilon 2
    dataframe['num_clusters'] = pd.Series(1, index=dataframe.index)
    evts = [i for i in range(0,len(dataframe))] #numbers events by dataframe index
    dataframe['event_number']=evts #renumbers events so they match dataframe indices
    events = dataframe.loc[dataframe['sum_tot']>800].index.to_numpy() #selection for events we separate tracks from
    col = dataframe['column']
    row = dataframe['row']

    for i in range(0,len(events)):
        data = np.array((list(col[events[i]]),list(row[events[i]])))
        data = data.T
        clusters = cluster.DBSCAN(eps=2).fit_predict(data) #algorithm outputs a list with an index corresponding to its cluster number
        df_clusters = pd.DataFrame(data=clusters,columns=['clust_val']) #make separate dataframe to access variables easily
        clust_indices = [df_clusters.loc[df_clusters['clust_val'] == np.unique(clusters)[i]].index.to_numpy() for i in range(0,len(np.unique(clusters)))] #entries are lists of indices corresponding to a given cluster
        for j in range(0,len(clust_indices)): #loop adds clusters as individual events to dataframe
            new_hit = pd.DataFrame({"event_number": events[i]+j, "timestamp": dataframe['timestamp'][events[i]], "npoints": len(clust_indices[j]), "column": [dataframe['column'][events[i]][clust_indices[j]]], "row": [dataframe['row'][events[i]][clust_indices[j]]], "BCID": [dataframe['BCID'][events[i]][clust_indices[j]]], "tot": [dataframe['tot'][events[i]][clust_indices[j]]], "sum_tot": np.sum(dataframe['tot'][events[i]][clust_indices[j]]), "num_clusters": len(clust_indices)},  index=[events[i]+j]) #adds new row for separated event from cluster
            dataframe = dataframe.append(new_hit, ignore_index=True,  sort=True)

    for i in range(0,len(events)): #drops the original multi-track events that were separated into clusters
        dataframe = dataframe.drop(dataframe.loc[dataframe.index==events[i]].index[0]) #removes multiple track events bypassing reassigned indices

    dataframe = dataframe.sort_values(by = ['timestamp']) #sorts updated from by event timestamp in ascending order
    dataframe.index = [i for i in range(0,len(dataframe))] #reindexes values from 0 to length of dataframe
    dataframe['event_number'] = dataframe.index #matches event numbers with indices

    return dataframe

def extract_parameters_from_h5_file(input_file):
    df_hits = pd.read_hdf(input_file, "Hits") #generates pandas dataframe of Hit-level data
    df_meta = pd.read_hdf(input_file, "meta_data") #generates dataframe of metadata
    df = pd.DataFrame()
    df['event_number'] = df_hits['event_number'].unique()
    events = df['event_number'].to_numpy()
    n = len(events)
    df['column'], df['row'], df['tot'], df['BCID'] = [df_hits['column'].loc[df_hits['event_number']==j].to_numpy() for j in events], [df_hits['row'].loc[df_hits['event_number']==j].to_numpy() for j in events], [df_hits['tot'].loc[df_hits['event_number']==j].to_numpy() for j in events], [df_hits['relative_BCID'].loc[df_hits['event_number']==j].to_numpy() for j in events]
    ts_entries = [df_meta['event_number'].loc[np.abs(df_meta['event_number']-events[j])<100] for j in range(0,n)]
    df['timestamp'] = [df_meta['timestamp_start'][ts_entries[j].index[0]] for j in range(0,n)]
    df['npoints'] = [len(df['column'][j]) for j in range(0,n)]
    df['sum_tot'] = [np.sum(df['tot'][j]) for j in range(0,n)]
    return df
'''
def quick_h5_extract(input_file):
    df_hits = pd.read_hdf(input_file, "Hits") #generates pandas dataframe of Hit-level data
    df_meta = pd.read_hdf(input_file, "meta_data") #generates dataframe of metadata
    df_meta['timestamp'] = df_meta['timestamp_start']
    df_meta = df_meta.drop(columns = ['error_code','timestamp_stop','timestamp_start'])
    df = df_meta.loc[df_meta['event_number'].isin(df_hits['event_number'])==True] #find events numbers that match
    df = df.loc[df.duplicated(['event_number']) == False] #remove duplicate numbers
    df.index = [i for i in range(0,len(df))]
    df_hits = df_hits.loc[df_hits['event_number'].isin(df_meta['event_number'])==True]
    df_hits.index = [i for i in range(0,len(df_hits))]
    events = df['event_number'].to_numpy()
    n = len(events)
    df['column'], df['row'], df['tot'], df['BCID'] = np.array(df_hits.groupby(['event_number'])['column'].apply(list)),np.array(df_hits.groupby(['event_number'])['row'].apply(list)),np.array(df_hits.groupby(['event_number'])['tot'].apply(list)),np.array(df_hits.groupby(['event_number'])['relative_BCID'].apply(list))
    df['column']=df['column'].apply(lambda x: np.array(x,dtype=np.uint8))
    df['row']=df['row'].apply(lambda x: np.array(x,dtype=np.uint16))
    df['tot']=df['tot'].apply(lambda x: np.array(x,dtype=np.uint8))
    df['BCID']=df['BCID'].apply(lambda x: np.array(x,dtype=np.uint8))
    df['npoints'] = [len(df['column'][j]) for j in range(0,n)]
    df = df.loc[df['npoints']<10000]
    df.index = [i for i in range(0,len(df))]
    df['sum_tot'] = [np.sum(df['tot'][j]) for j in range(0,len(df))]
    return df
'''
###Following scales fastest for large datasets
def quick_h5_extract(input_file):
    df_hits = pd.read_hdf(input_file, "Hits") #generates pandas dataframe of Hit-level data
    df_meta = pd.read_hdf(input_file, "meta_data") #generates dataframe of metadata
    df_meta['timestamp'] = df_meta['timestamp_start']
    df_meta = df_meta.drop(columns = ['error_code','timestamp_stop','timestamp_start'])
    df = df_meta.loc[df_meta['event_number'].isin(df_hits['event_number'])==True] #find events numbers that match
    df = df.loc[df.duplicated(['event_number']) == False] #remove duplicate numbers
    df.index = [i for i in range(0,len(df))]
    df_hits = df_hits.loc[df_hits['event_number'].isin(df_meta['event_number'])==True]
    df_hits.index = [i for i in range(0,len(df_hits))]
    events = df['event_number'].to_numpy()
    n = len(events)
    df_column = pd.DataFrame({'event_number': df_hits['event_number'], 'column': df_hits['column']})
    df_row = pd.DataFrame({'event_number': df_hits['event_number'], 'row': df_hits['row']})
    df_tot = pd.DataFrame({'event_number': df_hits['event_number'], 'tot': df_hits['tot']})
    df_BCID = pd.DataFrame({'event_number': df_hits['event_number'], 'BCID': df_hits['relative_BCID']})
    df['column'] = make_hits_arrays(df_column,'column')['column']
    df['column']=df['column'].apply(lambda x: np.array(x,dtype=np.uint8))
    df['row'] = make_hits_arrays(df_row,'row')['row']
    df['row']=df['row'].apply(lambda x: np.array(x,dtype=np.uint16))
    df['tot'] = make_hits_arrays(df_tot,'tot')['tot']
    df['tot']=df['tot'].apply(lambda x: np.array(x,dtype=np.uint8))
    df['BCID'] = make_hits_arrays(df_BCID,'BCID')['BCID']
    df['BCID']=df['BCID'].apply(lambda x: np.array(x,dtype=np.uint8))
    df['npoints'] = [len(df['column'][j]) for j in range(0,n)]
    df = df.loc[df['npoints']<10000]
    df.index = [i for i in range(0,len(df))]
    df['sum_tot'] = [np.sum(df['tot'][j]) for j in range(0,len(df))]
    return df

def make_hits_arrays(df,col):
         keys, values = df.sort_values('event_number').values.T
         ukeys, index = np.unique(keys, True)
         arrays = np.split(values, index[1:])
         df2 = pd.DataFrame({'event_number':ukeys, '%s'%(col):[np.array(a) for a in arrays]})
         return df2
