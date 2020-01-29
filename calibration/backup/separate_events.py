from root_pandas import read_root
import pandas as pd
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import sklearn.cluster as cluster
import array

def create_ntuple(f_input, f_output):
    output = ROOT.TFile(f_output, 'recreate')
    data = ROOT.TTree('data', 'data')
    input_df = read_dataframe(f_input) #dataframe from input root file
    df = separate_events(input_df) #dataframe with multi-track events separated into single events. Call df for short
    ### Variable names for ROOT branches###
    root_evt = array.array('i',[0])
    root_ts = array.array('d',[0])
    root_npoints = array.array('i',[0])
    root_sum_tot = array.array('i',[0])
    root_num_clusters = array.array('i',[0])
    root_col = array.array('i' , [0 for i in range(0,100000)])
    root_row = array.array('i' , [0 for i in range(0,100000)])
    root_tot = array.array('i' , [0 for i in range(0,100000)])
    root_BCID = array.array('i' , [0 for i in range(0,100000)])
    ### Define branches
    data.Branch("event_number", root_evt, "event_number/I")
    data.Branch("timestamp", root_ts , "timestamp/D")
    data.Branch("npoints", root_npoints , "npoints/I")
    data.Branch("column", root_col, "column[npoints]/I")
    data.Branch("row", root_row, "row[npoints]/I")
    data.Branch("BCID", root_BCID, "BCID[npoints]/I")
    data.Branch("tot", root_tot, "tot[npoints]/I")
    data.Branch("sum_tot", root_sum_tot, "sum_tot/I")
    data.Branch("num_clusters", root_num_clusters, "num_clusters/I")
    ### Variables from dataframe to be input into branches
    col, row, tot, bcid = df['column'], df['row'], df['tot'], df['BCID']
    timestamp = df['timestamp']
    events = df['event_number']
    npoints = df['npoints']
    tot_sum = df['sum_tot']
    num_clust = df['num_clusters']
    n = len(df)

    for j in range(0,n):
        root_evt[0] = events[j]
        root_ts[0] = timestamp[j]
        root_npoints[0] = npoints[j]
        root_sum_tot[0] = tot_sum[j]
        root_num_clusters[0] = num_clust[j]
        for i in range(0,root_npoints[0]):
            root_col[i] = col[j][i]-1
            root_row[i] = row[j][i]-1
            root_tot[i] = tot[j][i]
            root_BCID[i] = bcid[j][i]
        data.Fill()
    data.Print()
    output.Write()
    output.Close()

def separate_events(dataframe): #takes read in dataframe and separates multitrack events using DBSCAN with epsilon 2
    dataframe['num_clusters'] = pd.Series(1, index=dataframe.index)
    evts = [i for i in range(0,len(dataframe))] #numbers events by dataframe index
    dataframe['event_number']=evts #renumbers events so they match dataframe indices
    events = dataframe.loc[dataframe['sum_tot']>800].index #selection for events we separate tracks from
    col = dataframe['column']
    row = dataframe['row']

    for i in range(0,len(events)):
        data = np.array((list(col[events[i]]),list(row[events[i]])))
        data = data.T
        clusters = cluster.DBSCAN(eps=2).fit_predict(data) #algorithm outputs a list with an index corresponding to its cluster number
        df_clusters = pd.DataFrame(data=clusters,columns=['clust_val']) #make separate dataframe to access variables easily
        clust_indices = [df_clusters.loc[df_clusters['clust_val'] == np.unique(clusters)[i]].index for i in range(0,len(np.unique(clusters)))] #entries are lists of indices corresponding to a given cluster
        for j in range(0,len(clust_indices)): #loop adds clusters as individual events to dataframe
            new_hit = pd.DataFrame({"event_number": events[i]+j, "timestamp": dataframe['timestamp'][events[i]], "npoints": len(clust_indices[j]), "column": [dataframe['column'][events[i]][clust_indices[j]]], "row": [dataframe['row'][events[i]][clust_indices[j]]], "BCID": [dataframe['BCID'][events[i]][clust_indices[j]]], "tot": [dataframe['tot'][events[i]][clust_indices[j]]], "sum_tot": np.sum(dataframe['tot'][events[i]][clust_indices[j]]), "num_clusters": len(clust_indices)},  index=[events[i]+j]) #adds new row for separated event from cluster
            dataframe = dataframe.append(new_hit, ignore_index=True)

    for i in range(0,len(events)): #drops the original multi-track events that were separated into clusters
        dataframe = dataframe.drop(dataframe.loc[dataframe.index==events[i]].index[0]) #removes multiple track events bypassing reassigned indices

    dataframe = dataframe.sort_values(by = ['timestamp']) #sorts updated from by event timestamp in ascending order
    dataframe.index = [i for i in range(0,len(dataframe))] #reindexes values from 0 to length of dataframe
    dataframe['event_number'] = dataframe.index #matches event numbers with indeices

    return dataframe

def read_dataframe(f_input):
    
    dataframe = read_root(f_input) #input root file
    return dataframe
    

