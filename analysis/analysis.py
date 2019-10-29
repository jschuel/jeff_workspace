from root_pandas import read_root
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from array import array
import ROOT
import math
from os import sys

def make_multigraph(df):
    gr = make_plot(df, 789, 0)
    gr.SetMarkerStyle(20)
    gr.SetMarkerColor(4)
    
    gr1 = make_plot(df, 789, 1)
    gr1.SetMarkerStyle(20)
    gr1.SetMarkerColor(2)

    gr2 = make_plot(df, 789, -1)
    gr2.SetMarkerStyle(20)
    gr2.SetMarkerColor(2)
    
    gr3 = make_plot(df, 1576, 0)
    gr3.SetMarkerStyle(22)
    gr3.SetMarkerColor(64)
    
    gr4 = make_plot(df, 1576, 1)
    gr4.SetMarkerStyle(22)
    gr4.SetMarkerColor(96)

    gr5 = make_plot(df, 1576, -1)
    gr5.SetMarkerStyle(22)
    gr5.SetMarkerColor(96)
    
    gr6 = make_plot(df, 789, 2)
    gr6.SetMarkerStyle(20)
    gr6.SetMarkerColor(1)

    gr7 = make_plot(df, 789, -2)
    gr7.SetMarkerStyle(20)
    gr7.SetMarkerColor(1)

    gr8 = make_plot(df, 1576, 2)
    gr8.SetMarkerStyle(22)
    gr8.SetMarkerColor(1)

    gr9 = make_plot(df, 1576, -2)
    gr9.SetMarkerStyle(22)
    gr9.SetMarkerColor(1)

    mg = ROOT.TMultiGraph()
    mg.Add(gr)
    mg.Add(gr1)
    mg.Add(gr2)
    mg.Add(gr3)
    mg.Add(gr4)
    mg.Add(gr5)
    #mg.Add(gr6)
    #mg.Add(gr7)
    #mg.Add(gr8)
    #mg.Add(gr9)
    mg.Fit("pol1", "FQ")

    return mg

def make_plot(df, Nb, knob):
    index = get_indices(df, Nb, knob)
    y_plot = array('d' , df['y'][index])
    y_plot_err = array('d' , df['y_err'][index])
    x_plot = array('d', df['x'][index])
    x_plot_err = array('d', df['x_err'][index])
    if len(y_plot) > 0:
        gr = ROOT.TGraphErrors(len(y_plot), x_plot, y_plot, x_plot_err, y_plot_err)
        return gr
    else:
        gr = ROOT.TGraphErrors()
        return gr

def get_indices(dataframe, Nb, knob):
    #List indices where heuristic is not zero
    indices = dataframe.loc[(dataframe['y']>0) & (dataframe['knob']==knob) & (dataframe['Nb']==Nb)].index.to_numpy()
    return indices

def make_heuristic_dataframe(day, ring, module_id):
    df = read_root('/Users/vahsengrouplaptop/data/phase2/combined_SKB_TPC_ntuples/June_%s_%s.root'%(day, ring))

    bin_width = 120 #standard value
    nbins = math.floor(len(df)/bin_width)
    Z= 2.3

    #Create lists of numpy arrays of length bin_width. We will take the averages and st errors of these lists for analysis

    ts = []
    I = []
    P = []
    sigy = []
    Nb = []
    knob = []
    tpc_rate = []


    for ibin in range(0,nbins):
        ts.append(np.array([df['ts'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
        I.append(np.array([df['I'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
        P.append(np.array([df['P'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
        sigy.append(np.array([df['sigy'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
        Nb.append(np.array([df['Nb'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
        knob.append(np.array([df['knob'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
        tpc_rate.append(np.array([df['%s_neutrons'%(module_id)].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    
        #Compute averages and errors over bin_widths for each qty

    ts_avg = [np.mean(ts[i]) for i in range(0,len(ts))]
    ts_err = [np.std(ts[i])/math.sqrt(bin_width) for i in range(0,len(ts))]
    I_avg = [np.mean(I[i]) for i in range(0,len(I))]
    I_err = [np.std(I[i])/math.sqrt(bin_width) for i in range(0,len(I))]
    P_avg = [np.mean(P[i]) for i in range(0,len(P))]
    P_err = [np.std(P[i])/math.sqrt(bin_width) for i in range(0,len(P))]
    sigy_avg = [np.mean(sigy[i]) for i in range(0,len(sigy))]
    sigy_err = [np.std(sigy[i])/math.sqrt(bin_width) for i in range(0,len(sigy))]
    Nb_avg = [np.mean(Nb[i]) for i in range(0,len(Nb))]
    Nb_err = [np.std(Nb[i])/math.sqrt(bin_width) for i in range(0,len(Nb))]
    knob_avg = [np.mean(knob[i]) for i in range(0,len(knob))]
    knob_err = [np.std(knob[i])/math.sqrt(bin_width) for i in range(0,len(knob))]
    rate_avg = [np.mean(tpc_rate[i]) for i in range(0,len(tpc_rate))]
    rate_err = [np.std(tpc_rate[i])/math.sqrt(bin_width) for i in range(0,len(tpc_rate))]


    y = [rate_avg[i]/(I_avg[i]*P_avg[i]*Z**2) for i in range(0,len(ts_avg))]
    y_error = [y[i]*math.sqrt((rate_err[i]/rate_avg[i])**2 + (I_err[i]/I_avg[i])**2 + (P_err[i]/P_avg[i])**2) for i in range(0,len(ts_avg))]

    x = [I_avg[i]/(P_avg[i]*sigy_avg[i]*Nb_avg[i]*Z**2) for i in range(0,len(ts_avg))]
    x_err = [x[i]*math.sqrt((I_err[i]/I_avg[i])**2 + (P_err[i]/P_avg[i])**2 + (sigy_err[i]/sigy_avg[i])**2) for i in range (0,len(ts_avg))]

    #Create a dataframe so we can omit values whose averages are zero
    heuristic = pd.DataFrame({"x": x, "x_err": x_err, "y": y, "y_err": y_error, "Nb": Nb_avg, "knob": knob_avg})
    if ring == "HER":
        heuristic = heuristic.loc[((heuristic['Nb']==789) | (heuristic['Nb'] == 1576)) & ((heuristic['knob'] == 0) | (heuristic['knob'] == 1.0) | (heuristic['knob'] == 2.0) | (heuristic['knob'] == -1.0) | (heuristic['knob'] == -2.0))]
    if ring == "LER":
        heuristic = heuristic.loc[((heuristic['Nb']==789) | (heuristic['Nb'] == 1576)) & ((heuristic['knob'] == 0) | (heuristic['knob'] == 1.0) | (heuristic['knob'] == 1.9) | (heuristic['knob'] == -1.0) | (heuristic['knob'] == -1.9))]
    heuristic.index = [i for i in range(0,len(heuristic))]
    return heuristic

 
