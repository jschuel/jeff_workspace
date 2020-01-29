#cimport numpy as np
import numpy as np
import pandas as pd
import array
#from cpython cimport array
import ROOT

def make_ntuple(input_file, output_file):
    output = ROOT.TFile(output_file, 'recreate')
    global data
    data = ROOT.TTree('data', 'data')

    df_hits = pd.read_hdf(input_file, "Hits")
    df_meta = pd.read_hdf(input_file, "meta_data")
    events = df_hits['event_number'].unique()
    n = len(events)
    root_evt = array.array('i',[0])
    root_ts = array.array('d',[0])
    root_npoints = array.array('i',[0])
    root_sum_tot = array.array('i',[0])
    root_col = array.array('i' , [0 for i in range(0,100000)])
    root_row = array.array('i' , [0 for i in range(0,100000)])
    root_tot = array.array('i' , [0 for i in range(0,100000)])
    root_BCID = array.array('i' , [0 for i in range(0,100000)])
    data.Branch("event_number", root_evt, "event_number/I")
    data.Branch("timestamp", root_ts , "timestamp/D")
    data.Branch("npoints", root_npoints , "npoints/I")
    data.Branch("column", root_col, "column[npoints]/I")
    data.Branch("row", root_row, "row[npoints]/I")
    data.Branch("BCID", root_BCID, "BCID[npoints]/I")
    data.Branch("tot", root_tot, "tot[npoints]/I")
    data.Branch("sum_tot", root_sum_tot, "sum_tot/I")
    col, row, tot, bcid = [df_hits['column'].loc[df_hits['event_number']==j].to_numpy() for j in events], [df_hits['row'].loc[df_hits['event_number']==j].to_numpy() for j in events], [df_hits['tot'].loc[df_hits['event_number']==j].to_numpy() for j in events], [df_hits['relative_BCID'].loc[df_hits['event_number']==j].to_numpy() for j in events]
    ts_entries = [df_meta['event_number'].loc[np.abs(df_meta['event_number']-events[j])<100] for j in range(0,n)]
    timestamp = [df_meta['timestamp_start'][ts_entries[j].index[0]] for j in range(0,n)]
    npoints = [len(col[j]) for j in range(0,n)]
    tot_sum = [np.sum(tot[j]) for j in range(0,n)]

    for j in range(0,n):
        root_evt[0] = events[j]
        root_ts[0] = timestamp[j]
        root_npoints[0] = npoints[j]
        root_sum_tot[0] = tot_sum[j]
        for i in range(0,root_npoints[0]):
            root_col[i] = col[j][i]-1
            root_row[i] = row[j][i]-1
            root_tot[i] = tot[j][i]
            root_BCID[i] = bcid[j][i]
        data.Fill()

    data.Print()
    output.Write()
    output.Close()
