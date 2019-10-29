from root_pandas import read_root
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from array import array
import ROOT
import math
from os import sys

day = sys.argv[1]
ring = sys.argv[2]
module_id = sys.argv[3]

df = read_root('/Users/vahsengrouplaptop/data/phase2/combined_SKB_TPC_ntuples/June_%s_%s.root'%(day, ring))

bin_width = 300
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
x_err = [x[i]*math.sqrt((I_err[i]/I_avg[i])**2 + (P_err[i]/P_avg[i])**2 + (sigy_err[i]/sigy_avg[i])**2) for i in range (0,len(ts))]

#Create a dataframe so we can omit values whose averages are zero
heuristic = pd.DataFrame({"x": x, "x_err": x_err, "y": y, "y_err": y_error, "Nb": Nb_avg, "knob": knob_avg})
if ring == "HER":
    heuristic = heuristic.loc[((heuristic['Nb']==789) | (heuristic['Nb'] == 1576)) & ((heuristic['knob'] == 0) | (heuristic['knob'] == 1.0) | (heuristic['knob'] == 2.0) | (heuristic['knob'] == -1.0) | (heuristic['knob'] == -2.0))]
if ring == "LER":
    heuristic = heuristic.loc[((heuristic['Nb']==789) | (heuristic['Nb'] == 1576)) & ((heuristic['knob'] == 0) | (heuristic['knob'] == 1.0) | (heuristic['knob'] == 1.9) | (heuristic['knob'] == -1.0) | (heuristic['knob'] == -1.9))]
heuristic.index = [i for i in range(0,len(heuristic))]

#List indices where heuristic is not zero
indices = heuristic.loc[(heuristic['y']>0)].index.to_numpy()

#indices789 = heuristic.loc[(heuristic['y']>0) & (heuristic['Nb']==789)].index.to_numpy()

#indices1576 = heuristic.loc[(heuristic['y']>0) & (heuristic['Nb']==1576)].index.to_numpy()

indices_knob0_789 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']==0) & (heuristic['Nb']==789)].index.to_numpy()

indices_knob1_789 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']==1) & (heuristic['Nb']==789)].index.to_numpy()

indices_knobneg1_789 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']==-1) & (heuristic['Nb']==789)].index.to_numpy()

indices_knob2_789 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']>1) & (heuristic['Nb']==789)].index.to_numpy()

indices_knobneg2_789 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']<-1) & (heuristic['Nb']==789)].index.to_numpy()

indices_knob0_1576 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']==0) & (heuristic['Nb']==1576)].index.to_numpy()

indices_knob1_1576 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']==1) & (heuristic['Nb']==1576)].index.to_numpy()

indices_knobneg1_1576 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']==-1) & (heuristic['Nb']==1576)].index.to_numpy()

indices_knob2_1576 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']>1) & (heuristic['Nb']==1576)].index.to_numpy()

indices_knobneg2_1576 = heuristic.loc[(heuristic['y']>0) & (heuristic['knob']<-1) & (heuristic['Nb']==1576)].index.to_numpy()

### DEFINE VARIABLES FOR ROOT PLOTS###

y_plot = array('d' , heuristic['y'][indices])
y_plot_err = array('d' , heuristic['y_err'][indices])
x_plot = array('d', heuristic['x'][indices])
x_plot_err = array('d', heuristic['x_err'][indices])

y_plot_knob0_789 = array('d' , heuristic['y'][indices_knob0_789])
y_plot_err_knob0_789 = array('d' , heuristic['y_err'][indices_knob0_789])
x_plot_knob0_789 = array('d', heuristic['x'][indices_knob0_789])
x_plot_err_knob0_789 = array('d', heuristic['x_err'][indices_knob0_789])

y_plot_knob1_789 = array('d' , heuristic['y'][indices_knob1_789])
y_plot_err_knob1_789 = array('d' , heuristic['y_err'][indices_knob1_789])
x_plot_knob1_789 = array('d', heuristic['x'][indices_knob1_789])
x_plot_err_knob1_789 = array('d', heuristic['x_err'][indices_knob1_789])

y_plot_knob2_789 = array('d' , heuristic['y'][indices_knob2_789])
y_plot_err_knob2_789 = array('d' , heuristic['y_err'][indices_knob2_789])
x_plot_knob2_789 = array('d', heuristic['x'][indices_knob2_789])
x_plot_err_knob2_789 = array('d', heuristic['x_err'][indices_knob2_789])

y_plot_knobneg1_789 = array('d' , heuristic['y'][indices_knobneg1_789])
y_plot_err_knobneg1_789 = array('d' , heuristic['y_err'][indices_knobneg1_789])
x_plot_knobneg1_789 = array('d', heuristic['x'][indices_knobneg1_789])
x_plot_err_knobneg1_789 = array('d', heuristic['x_err'][indices_knobneg1_789])

y_plot_knobneg2_789 = array('d' , heuristic['y'][indices_knobneg2_789])
y_plot_err_knobneg2_789 = array('d' , heuristic['y_err'][indices_knobneg2_789])
x_plot_knobneg2_789 = array('d', heuristic['x'][indices_knobneg2_789])
x_plot_err_knobneg2_789 = array('d', heuristic['x_err'][indices_knobneg2_789])

y_plot_knob0_1576 = array('d' , heuristic['y'][indices_knob0_1576])
y_plot_err_knob0_1576 = array('d' , heuristic['y_err'][indices_knob0_1576])
x_plot_knob0_1576 = array('d', heuristic['x'][indices_knob0_1576])
x_plot_err_knob0_1576 = array('d', heuristic['x_err'][indices_knob0_1576])

y_plot_knob1_1576 = array('d' , heuristic['y'][indices_knob1_1576])
y_plot_err_knob1_1576 = array('d' , heuristic['y_err'][indices_knob1_1576])
x_plot_knob1_1576 = array('d', heuristic['x'][indices_knob1_1576])
x_plot_err_knob1_1576 = array('d', heuristic['x_err'][indices_knob1_1576])

y_plot_knob2_1576 = array('d' , heuristic['y'][indices_knob2_1576])
y_plot_err_knob2_1576 = array('d' , heuristic['y_err'][indices_knob2_1576])
x_plot_knob2_1576 = array('d', heuristic['x'][indices_knob2_1576])
x_plot_err_knob2_1576 = array('d', heuristic['x_err'][indices_knob2_1576])

y_plot_knobneg1_1576 = array('d' , heuristic['y'][indices_knobneg1_1576])
y_plot_err_knobneg1_1576 = array('d' , heuristic['y_err'][indices_knobneg1_1576])
x_plot_knobneg1_1576 = array('d', heuristic['x'][indices_knobneg1_1576])
x_plot_err_knobneg1_1576 = array('d', heuristic['x_err'][indices_knobneg1_1576])

y_plot_knobneg2_1576 = array('d' , heuristic['y'][indices_knobneg2_1576])
y_plot_err_knobneg2_1576 = array('d' , heuristic['y_err'][indices_knobneg2_1576])
x_plot_knobneg2_1576 = array('d', heuristic['x'][indices_knobneg2_1576])
x_plot_err_knobneg2_1576 = array('d', heuristic['x_err'][indices_knobneg2_1576])

##Plot ROOT plots##

c1 = ROOT.TCanvas("c1", "heuristic fits", 1024, 768)
c1.Divide(1,2)
c1.cd(1)

f_lin = ROOT.TF1("f_lin", "[0]+[1]*x" )

gr_tpc = ROOT.TGraphErrors(len(indices), x_plot, y_plot, x_plot_err, y_plot_err)
gr_tpc.SetMarkerColor(1)
gr_tpc.SetMarkerStyle(20)
gr_tpc.SetMarkerSize(1)
fit_tpc = gr_tpc.Fit("f_lin", "S")
#gr_tpc.GetFunction("f_lin").SetLineColor(2)

gr_tpc_knob0_789 = ROOT.TGraphErrors(len(indices_knob0_789), x_plot_knob0_789, y_plot_knob0_789, x_plot_err_knob0_789, y_plot_err_knob0_789)
gr_tpc_knob0_789.SetMarkerColor(4)
gr_tpc_knob0_789.SetMarkerStyle(20)
gr_tpc_knob0_789.SetMarkerSize(1)

gr_tpc_knob1_789 = ROOT.TGraphErrors(len(indices_knob1_789), x_plot_knob1_789, y_plot_knob1_789, x_plot_err_knob1_789, y_plot_err_knob1_789)
gr_tpc_knob1_789.SetMarkerColor(2)
gr_tpc_knob1_789.SetMarkerStyle(20)
gr_tpc_knob1_789.SetMarkerSize(1)

gr_tpc_knob2_789 = ROOT.TGraphErrors(len(indices_knob2_789), x_plot_knob2_789, y_plot_knob2_789, x_plot_err_knob2_789, y_plot_err_knob2_789)
gr_tpc_knob2_789.SetMarkerColor(1)
gr_tpc_knob2_789.SetMarkerStyle(20)
gr_tpc_knob2_789.SetMarkerSize(1)

gr_tpc_knobneg1_789 = ROOT.TGraphErrors(len(indices_knobneg1_789), x_plot_knobneg1_789, y_plot_knobneg1_789, x_plot_err_knobneg1_789, y_plot_err_knobneg1_789)
gr_tpc_knobneg1_789.SetMarkerColor(2)
gr_tpc_knobneg1_789.SetMarkerStyle(20)
gr_tpc_knobneg1_789.SetMarkerSize(1)

gr_tpc_knobneg2_789 = ROOT.TGraphErrors(len(indices_knobneg2_789), x_plot_knobneg2_789, y_plot_knobneg2_789, x_plot_err_knobneg2_789, y_plot_err_knobneg2_789)
gr_tpc_knobneg2_789.SetMarkerColor(1)
gr_tpc_knobneg2_789.SetMarkerStyle(20)
gr_tpc_knobneg2_789.SetMarkerSize(1)

gr_tpc_knob0_1576 = ROOT.TGraphErrors(len(indices_knob0_1576), x_plot_knob0_1576, y_plot_knob0_1576, x_plot_err_knob0_1576, y_plot_err_knob0_1576)
gr_tpc_knob0_1576.SetMarkerColor(4)
gr_tpc_knob0_1576.SetMarkerStyle(21)
gr_tpc_knob0_1576.SetMarkerSize(1)

gr_tpc_knob1_1576 = ROOT.TGraphErrors(len(indices_knob1_1576), x_plot_knob1_1576, y_plot_knob1_1576, x_plot_err_knob1_1576, y_plot_err_knob1_1576)
gr_tpc_knob1_1576.SetMarkerColor(2)
gr_tpc_knob1_1576.SetMarkerStyle(21)
gr_tpc_knob1_1576.SetMarkerSize(1)

gr_tpc_knob2_1576 = ROOT.TGraphErrors(len(indices_knob2_1576), x_plot_knob2_1576, y_plot_knob2_1576, x_plot_err_knob2_1576, y_plot_err_knob2_1576)
gr_tpc_knob2_1576.SetMarkerColor(1)
gr_tpc_knob2_1576.SetMarkerStyle(21)
gr_tpc_knob2_1576.SetMarkerSize(1)
'''
gr_tpc_knobneg1_1576 = ROOT.TGraphErrors(len(indices_knobneg1_1576), x_plot_knobneg1_1576, y_plot_knobneg1_1576, x_plot_err_knobneg1_1576, y_plot_err_knobneg1_1576)
gr_tpc_knobneg1_1576.SetMarkerColor(2)
gr_tpc_knobneg1_1576.SetMarkerStyle(21)
gr_tpc_knobneg1_1576.SetMarkerSize(1)

gr_tpc_knobneg2_1576 = ROOT.TGraphErrors(len(indices_knobneg2_1576), x_plot_knobneg2_1576, y_plot_knobneg2_1576, x_plot_err_knobneg2_1576, y_plot_err_knobneg2_1576)
gr_tpc_knobneg2_1576.SetMarkerColor(1)
gr_tpc_knobneg2_1576.SetMarkerStyle(21)
gr_tpc_knobneg2_1576.SetMarkerSize(1)
'''

mg_BWD = ROOT.TMultiGraph()
mg_BWD.Add(gr_tpc)
mg_BWD.Add(gr_tpc_knob0_789)
mg_BWD.Add(gr_tpc_knob1_789)
mg_BWD.Add(gr_tpc_knob2_789)
mg_BWD.Add(gr_tpc_knobneg1_789)
mg_BWD.Add(gr_tpc_knobneg2_789)
mg_BWD.Add(gr_tpc_knob0_1576)
mg_BWD.Add(gr_tpc_knob1_1576)
mg_BWD.Add(gr_tpc_knob2_1576)
#mg_BWD.Add(gr_tpc_knobneg1_1576)
#mg_BWD.Add(gr_tpc_knobneg2_1576)
#mg_BWD.GetXaxis().SetLimits(0,100000)
mg_BWD.Draw("AP")

l_BWD = ROOT.TLegend(0.1,0.7,0.28,0.6)
l_BWD.AddEntry(gr_tpc, "%s"%(module_id), "Ep")
l_BWD.Draw()

###Make histogram of rates using B and T coefficients

B_HER_tpc = 3.4e2
B_LER_tpc = 8.7e1
T_HER_tpc = 0
T_LER_tpc = 4.4e-5


#HER_beam_gas_elepaio = B_HER_elepaio * 320 * 1.6e-8 * Z**2
#LER_beam_gas_elepaio = B_LER_elepaio * 350 * 6.4e-8 * Z**2
#HER_touschek_elepaio = T_HER_elepaio * 320**2/(60*789)
#LER_touschek_elepaio = T_LER_elepaio * 350**2/(130*789)

#print("Total rates for elepaio are, %s, %s, %s, %s"%(HER_beam_gas_elepaio, LER_beam_gas_elepaio, HER_touschek_elepaio, LER_touschek_elepaio))  
