from root_pandas import read_root
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from array import array
import ROOT
import math
from os import sys
from compute_heuristic_y import *

day = sys.argv[1]
ring = sys.argv[2]

df = read_root('/Users/vahsengrouplaptop/data/phase3/combined_SKB_TPC_ntuples/May_%s_%s.root'%(day, ring))

bin_width = 450
nbins = math.floor(len(df)/bin_width)
Z= 2.3

#Create lists of numpy arrays of length bin_width. We will take the averages and st errors of these lists for analysis

ts = []
I = []
P = []
sigy = []
Nb = []
elepaio_rate = []
palila_rate = []
tako_rate = []
humu_rate = []
nene_rate = []

for ibin in range(0,nbins):
    ts.append(np.array([df['ts'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    I.append(np.array([df['I'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    P.append(np.array([df['P'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    sigy.append(np.array([df['sigy'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    Nb.append(np.array([df['Nb'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    elepaio_rate.append(np.array([df['elepaio_neutrons'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    palila_rate.append(np.array([df['palila_neutrons'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    tako_rate.append(np.array([df['tako_neutrons'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    humu_rate.append(np.array([df['humu_neutrons'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    nene_rate.append(np.array([df['nene_neutrons'].to_numpy()[(ibin*bin_width)+i] for i in range(0,bin_width)]))
    
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
elepaio_rate_avg = [np.mean(elepaio_rate[i]) for i in range(0,len(elepaio_rate))]
elepaio_rate_err = [np.std(elepaio_rate[i])/math.sqrt(bin_width) for i in range(0,len(elepaio_rate))]
palila_rate_avg = [np.mean(palila_rate[i]) for i in range(0,len(palila_rate))]
palila_rate_err = [np.std(palila_rate[i])/math.sqrt(bin_width) for i in range(0,len(palila_rate))]
tako_rate_avg = [np.mean(tako_rate[i]) for i in range(0,len(tako_rate))]
tako_rate_err = [np.std(tako_rate[i])/math.sqrt(bin_width) for i in range(0,len(tako_rate))]
humu_rate_avg = [np.mean(humu_rate[i]) for i in range(0,len(humu_rate))]
humu_rate_err = [np.std(humu_rate[i])/math.sqrt(bin_width) for i in range(0,len(humu_rate))]
nene_rate_avg = [np.mean(nene_rate[i]) for i in range(0,len(nene_rate))]
nene_rate_err = [np.std(nene_rate[i])/math.sqrt(bin_width) for i in range(0,len(nene_rate))]


y_elepaio = [y_heuristic(elepaio_rate_avg[i], I_avg[i], P_avg[i], Z) for i in range(0,len(ts_avg))]
y_elepaio_error = [y_heuristic_error(y_elepaio[i], elepaio_rate_avg[i], elepaio_rate_err[i], I_avg[i], I_err[i], P_avg[i], P_err[i]) for i in range(0,len(ts_avg))]
y_palila = [y_heuristic(palila_rate_avg[i], I_avg[i], P_avg[i], Z) for i in range(0,len(ts_avg))]
y_palila_error = [y_heuristic_error(y_palila[i], palila_rate_avg[i], palila_rate_err[i], I_avg[i], I_err[i], P_avg[i], P_err[i]) for i in range(0,len(ts_avg))]
y_tako = [y_heuristic(tako_rate_avg[i], I_avg[i], P_avg[i], Z) for i in range(0,len(ts_avg))]
y_tako_error = [y_heuristic_error(y_tako[i], tako_rate_avg[i], tako_rate_err[i], I_avg[i], I_err[i], P_avg[i], P_err[i]) for i in range(0,len(ts_avg))]
y_humu = [y_heuristic(humu_rate_avg[i], I_avg[i], P_avg[i], Z) for i in range(0,len(ts_avg))]
y_humu_error = [y_heuristic_error(y_humu[i], humu_rate_avg[i], humu_rate_err[i], I_avg[i], I_err[i], P_avg[i], P_err[i]) for i in range(0,len(ts_avg))]
y_nene = [y_heuristic(nene_rate_avg[i], I_avg[i], P_avg[i], Z) for i in range(0,len(ts_avg))]
y_nene_error = [y_heuristic_error(y_nene[i], nene_rate_avg[i], nene_rate_err[i], I_avg[i], I_err[i], P_avg[i], P_err[i]) for i in range(0,len(ts_avg))]

x = [I_avg[i]/(P_avg[i]*sigy_avg[i]*Nb_avg[i]*Z**2) for i in range(0,len(ts_avg))]
x_err = [x[i]*math.sqrt((I_err[i]/I_avg[i])**2 + (P_err[i]/P_avg[i])**2 + (sigy_err[i]/sigy_avg[i])**2) for i in range (0,len(ts))]

#Create a dataframe so we can omit values whose averages are zero
heuristic = pd.DataFrame({"x": x, "x_err": x_err, "elepaio": y_elepaio, "elepaio_error": y_elepaio_error, "palila": y_palila, "palila_error": y_palila_error, "tako": y_tako, "tako_error": y_tako_error, "humu": y_humu, "humu_error": y_humu_error, "nene": y_nene, "nene_error": y_nene_error, "Nb": Nb_avg})
heuristic = heuristic.loc[(heuristic['Nb']==789) | (heuristic['Nb'] == 1576) | (heuristic['Nb'] == 395)]
heuristic.index = [i for i in range(0,len(heuristic))]

#List indices where heuristic is not zero
elepaio_indices = heuristic.loc[(heuristic['elepaio']>0)].index.to_numpy()
palila_indices = heuristic.loc[(heuristic['palila']>0)].index.to_numpy()
tako_indices = heuristic.loc[(heuristic['tako']>0)].index.to_numpy()
humu_indices = heuristic.loc[(heuristic['humu']>0)].index.to_numpy()
nene_indices = heuristic.loc[(heuristic['nene']>0)].index.to_numpy()

elepaio_indices395 = heuristic.loc[(heuristic['elepaio']>0) & (heuristic['Nb']==395)].index.to_numpy()
palila_indices395 = heuristic.loc[(heuristic['palila']>0) & (heuristic['Nb']==395)].index.to_numpy()
tako_indices395 = heuristic.loc[(heuristic['tako']>0) & (heuristic['Nb']==395)].index.to_numpy()
humu_indices395 = heuristic.loc[(heuristic['humu']>0) & (heuristic['Nb']==395)].index.to_numpy()
nene_indices395 = heuristic.loc[(heuristic['nene']>0) & (heuristic['Nb']==395)].index.to_numpy()

elepaio_indices789 = heuristic.loc[(heuristic['elepaio']>0) & (heuristic['Nb']==789)].index.to_numpy()
palila_indices789 = heuristic.loc[(heuristic['palila']>0) & (heuristic['Nb']==789)].index.to_numpy()
tako_indices789 = heuristic.loc[(heuristic['tako']>0) & (heuristic['Nb']==789)].index.to_numpy()
humu_indices789 = heuristic.loc[(heuristic['humu']>0) & (heuristic['Nb']==789)].index.to_numpy()
nene_indices789 = heuristic.loc[(heuristic['nene']>0) & (heuristic['Nb']==789)].index.to_numpy()

elepaio_indices1576 = heuristic.loc[(heuristic['elepaio']>0) & (heuristic['Nb']==1576)].index.to_numpy()
palila_indices1576 = heuristic.loc[(heuristic['palila']>0) & (heuristic['Nb']==1576)].index.to_numpy()
tako_indices1576 = heuristic.loc[(heuristic['tako']>0) & (heuristic['Nb']==1576)].index.to_numpy()
humu_indices1576 = heuristic.loc[(heuristic['humu']>0) & (heuristic['Nb']==1576)].index.to_numpy()
nene_indices1576 = heuristic.loc[(heuristic['nene']>0) & (heuristic['Nb']==1576)].index.to_numpy()
### DEFINE VARIABLES FOR ROOT PLOTS###

elepaio_y_plot = array('d' , heuristic['elepaio'][elepaio_indices])
elepaio_y_plot_err = array('d' , heuristic['elepaio_error'][elepaio_indices])
elepaio_x_plot = array('d', heuristic['x'][elepaio_indices])
elepaio_x_plot_err = array('d', heuristic['x_err'][elepaio_indices])

palila_y_plot = array('d' , heuristic['palila'][palila_indices])
palila_y_plot_err = array('d' , heuristic['palila_error'][palila_indices])
palila_x_plot = array('d', heuristic['x'][palila_indices])
palila_x_plot_err = array('d', heuristic['x_err'][palila_indices])

tako_y_plot = array('d' , heuristic['tako'][tako_indices])
tako_y_plot_err = array('d' , heuristic['tako_error'][tako_indices])
tako_x_plot = array('d', heuristic['x'][tako_indices])
tako_x_plot_err = array('d', heuristic['x_err'][tako_indices])

humu_y_plot = array('d' , heuristic['humu'][humu_indices])
humu_y_plot_err = array('d' , heuristic['humu_error'][humu_indices])
humu_x_plot = array('d', heuristic['x'][humu_indices])
humu_x_plot_err = array('d', heuristic['x_err'][humu_indices])

nene_y_plot = array('d' , heuristic['nene'][nene_indices])
nene_y_plot_err = array('d' , heuristic['nene_error'][nene_indices])
nene_x_plot = array('d', heuristic['x'][nene_indices])
nene_x_plot_err = array('d', heuristic['x_err'][nene_indices])

elepaio_y_plot395 = array('d' , heuristic['elepaio'][elepaio_indices395])
elepaio_y_plot_err395 = array('d' , heuristic['elepaio_error'][elepaio_indices395])
elepaio_x_plot395 = array('d', heuristic['x'][elepaio_indices395])
elepaio_x_plot_err395 = array('d', heuristic['x_err'][elepaio_indices395])

palila_y_plot395 = array('d' , heuristic['palila'][palila_indices395])
palila_y_plot_err395 = array('d' , heuristic['palila_error'][palila_indices395])
palila_x_plot395 = array('d', heuristic['x'][palila_indices395])
palila_x_plot_err395 = array('d', heuristic['x_err'][palila_indices395])

tako_y_plot395 = array('d' , heuristic['tako'][tako_indices395])
tako_y_plot_err395 = array('d' , heuristic['tako_error'][tako_indices395])
tako_x_plot395 = array('d', heuristic['x'][tako_indices395])
tako_x_plot_err395 = array('d', heuristic['x_err'][tako_indices395])

humu_y_plot395 = array('d' , heuristic['humu'][humu_indices395])
humu_y_plot_err395 = array('d' , heuristic['humu_error'][humu_indices395])
humu_x_plot395 = array('d', heuristic['x'][humu_indices395])
humu_x_plot_err395 = array('d', heuristic['x_err'][humu_indices395])

nene_y_plot395 = array('d' , heuristic['nene'][nene_indices395])
nene_y_plot_err395 = array('d' , heuristic['nene_error'][nene_indices395])
nene_x_plot395 = array('d', heuristic['x'][nene_indices395])
nene_x_plot_err395 = array('d', heuristic['x_err'][nene_indices395])

elepaio_y_plot789 = array('d' , heuristic['elepaio'][elepaio_indices789])
elepaio_y_plot_err789 = array('d' , heuristic['elepaio_error'][elepaio_indices789])
elepaio_x_plot789 = array('d', heuristic['x'][elepaio_indices789])
elepaio_x_plot_err789 = array('d', heuristic['x_err'][elepaio_indices789])

palila_y_plot789 = array('d' , heuristic['palila'][palila_indices789])
palila_y_plot_err789 = array('d' , heuristic['palila_error'][palila_indices789])
palila_x_plot789 = array('d', heuristic['x'][palila_indices789])
palila_x_plot_err789 = array('d', heuristic['x_err'][palila_indices789])

tako_y_plot789 = array('d' , heuristic['tako'][tako_indices789])
tako_y_plot_err789 = array('d' , heuristic['tako_error'][tako_indices789])
tako_x_plot789 = array('d', heuristic['x'][tako_indices789])
tako_x_plot_err789 = array('d', heuristic['x_err'][tako_indices789])

humu_y_plot789 = array('d' , heuristic['humu'][humu_indices789])
humu_y_plot_err789 = array('d' , heuristic['humu_error'][humu_indices789])
humu_x_plot789 = array('d', heuristic['x'][humu_indices789])
humu_x_plot_err789 = array('d', heuristic['x_err'][humu_indices789])

nene_y_plot789 = array('d' , heuristic['nene'][nene_indices789])
nene_y_plot_err789 = array('d' , heuristic['nene_error'][nene_indices789])
nene_x_plot789 = array('d', heuristic['x'][nene_indices789])
nene_x_plot_err789 = array('d', heuristic['x_err'][nene_indices789])

elepaio_y_plot1576 = array('d' , heuristic['elepaio'][elepaio_indices1576])
elepaio_y_plot_err1576 = array('d' , heuristic['elepaio_error'][elepaio_indices1576])
elepaio_x_plot1576 = array('d', heuristic['x'][elepaio_indices1576])
elepaio_x_plot_err1576 = array('d', heuristic['x_err'][elepaio_indices1576])

palila_y_plot1576 = array('d' , heuristic['palila'][palila_indices1576])
palila_y_plot_err1576 = array('d' , heuristic['palila_error'][palila_indices1576])
palila_x_plot1576 = array('d', heuristic['x'][palila_indices1576])
palila_x_plot_err1576 = array('d', heuristic['x_err'][palila_indices1576])

tako_y_plot1576 = array('d' , heuristic['tako'][tako_indices1576])
tako_y_plot_err1576 = array('d' , heuristic['tako_error'][tako_indices1576])
tako_x_plot1576 = array('d', heuristic['x'][tako_indices1576])
tako_x_plot_err1576 = array('d', heuristic['x_err'][tako_indices1576])

humu_y_plot1576 = array('d' , heuristic['humu'][humu_indices1576])
humu_y_plot_err1576 = array('d' , heuristic['humu_error'][humu_indices1576])
humu_x_plot1576 = array('d', heuristic['x'][humu_indices1576])
humu_x_plot_err1576 = array('d', heuristic['x_err'][humu_indices1576])

nene_y_plot1576 = array('d' , heuristic['nene'][nene_indices1576])
nene_y_plot_err1576 = array('d' , heuristic['nene_error'][nene_indices1576])
nene_x_plot1576 = array('d', heuristic['x'][nene_indices1576])
nene_x_plot_err1576 = array('d', heuristic['x_err'][nene_indices1576])

##Plot ROOT plots##

c1 = ROOT.TCanvas("c1", "heuristic fits", 1024, 768)
c1.Divide(1,2)
c1.cd(1)

f_lin = ROOT.TF1("f_lin", "[0]+[1]*x" )

gr_elepaio = ROOT.TGraphErrors(len(elepaio_indices), elepaio_x_plot, elepaio_y_plot, elepaio_x_plot_err, elepaio_y_plot_err)
gr_elepaio.SetMarkerColor(1)
gr_elepaio.SetMarkerStyle(24)
gr_elepaio.SetMarkerSize(1)
fit_elepaio = gr_elepaio.Fit("f_lin", "S")
#gr_elepaio.GetFunction("f_lin").SetLineColor(2)

gr_elepaio395 = ROOT.TGraphErrors(len(elepaio_indices395), elepaio_x_plot395, elepaio_y_plot395, elepaio_x_plot_err395, elepaio_y_plot_err395)
gr_elepaio395.SetMarkerColor(6)
gr_elepaio395.SetMarkerStyle(24)
gr_elepaio395.SetMarkerSize(1)

gr_elepaio789 = ROOT.TGraphErrors(len(elepaio_indices789), elepaio_x_plot789, elepaio_y_plot789, elepaio_x_plot_err789, elepaio_y_plot_err789)
gr_elepaio789.SetMarkerColor(4)
gr_elepaio789.SetMarkerStyle(24)
gr_elepaio789.SetMarkerSize(1)

gr_elepaio1576 = ROOT.TGraphErrors(len(elepaio_indices1576), elepaio_x_plot1576, elepaio_y_plot1576, elepaio_x_plot_err1576, elepaio_y_plot_err1576)
gr_elepaio1576.SetMarkerColor(8)
gr_elepaio1576.SetMarkerStyle(24)
gr_elepaio1576.SetMarkerSize(1)

gr_palila = ROOT.TGraphErrors(len(palila_indices), palila_x_plot, palila_y_plot, palila_x_plot_err, palila_y_plot_err)
gr_palila.SetMarkerColor(1)
gr_palila.SetMarkerStyle(21)
gr_palila.SetMarkerSize(1)
fit_palila = gr_palila.Fit("f_lin", "S")
#gr_palila.GetFunction("f_lin").SetLineColor(2)

gr_palila395 = ROOT.TGraphErrors(len(palila_indices395), palila_x_plot395, palila_y_plot395, palila_x_plot_err395, palila_y_plot_err395)
gr_palila395.SetMarkerColor(6)
gr_palila395.SetMarkerStyle(21)
gr_palila395.SetMarkerSize(1)

gr_palila789 = ROOT.TGraphErrors(len(palila_indices789), palila_x_plot789, palila_y_plot789, palila_x_plot_err789, palila_y_plot_err789)
gr_palila789.SetMarkerColor(4)
gr_palila789.SetMarkerStyle(21)
gr_palila789.SetMarkerSize(1)

gr_palila1576 = ROOT.TGraphErrors(len(palila_indices1576), palila_x_plot1576, palila_y_plot1576, palila_x_plot_err1576, palila_y_plot_err1576)
gr_palila1576.SetMarkerColor(8)
gr_palila1576.SetMarkerStyle(21)
gr_palila1576.SetMarkerSize(1)

gr_tako = ROOT.TGraphErrors(len(tako_indices), tako_x_plot, tako_y_plot, tako_x_plot_err, tako_y_plot_err)
gr_tako.SetMarkerColor(1)
gr_tako.SetMarkerStyle(22)
gr_tako.SetMarkerSize(1)
fit_tako = gr_tako.Fit("f_lin", "S")
#gr_tako.GetFunction("f_lin").SetLineColor(1)

gr_tako395 = ROOT.TGraphErrors(len(tako_indices395), tako_x_plot395, tako_y_plot395, tako_x_plot_err395, tako_y_plot_err395)
gr_tako395.SetMarkerColor(6)
gr_tako395.SetMarkerStyle(22)
gr_tako395.SetMarkerSize(1)

gr_tako789 = ROOT.TGraphErrors(len(tako_indices789), tako_x_plot789, tako_y_plot789, tako_x_plot_err789, tako_y_plot_err789)
gr_tako789.SetMarkerColor(4)
gr_tako789.SetMarkerStyle(22)
gr_tako789.SetMarkerSize(1)

gr_tako1576 = ROOT.TGraphErrors(len(tako_indices1576), tako_x_plot1576, tako_y_plot1576, tako_x_plot_err1576, tako_y_plot_err1576)
gr_tako1576.SetMarkerColor(8)
gr_tako1576.SetMarkerStyle(22)
gr_tako1576.SetMarkerSize(1)

mg_BWD = ROOT.TMultiGraph()
mg_BWD.Add(gr_elepaio)
mg_BWD.Add(gr_elepaio395)
mg_BWD.Add(gr_elepaio789)
mg_BWD.Add(gr_elepaio1576)
mg_BWD.Add(gr_palila)
mg_BWD.Add(gr_palila395)
mg_BWD.Add(gr_palila789)
mg_BWD.Add(gr_palila1576)
mg_BWD.Add(gr_tako)
mg_BWD.Add(gr_tako395)
mg_BWD.Add(gr_tako789)
mg_BWD.Add(gr_tako1576)
#mg_BWD.GetXaxis().SetLimits(0,100000)
mg_BWD.Draw("AP")

l_BWD = ROOT.TLegend(0.1,0.7,0.28,0.6)
l_BWD.AddEntry(gr_elepaio, "Elepaio", "Ep")
l_BWD.AddEntry(gr_palila, "Palila", "Ep")
l_BWD.AddEntry(gr_tako, "Tako", "Ep")
l_BWD.Draw()

c1.Update()

c1.cd(2)

gr_humu = ROOT.TGraphErrors(len(humu_indices), humu_x_plot, humu_y_plot, humu_x_plot_err, humu_y_plot_err)
gr_humu.SetMarkerColor(1)
gr_humu.SetMarkerStyle(21)
gr_humu.SetMarkerSize(1)
fit_humu = gr_humu.Fit("f_lin", "S")
#gr_humu.GetFunction("f_lin").SetLineColor(1)

gr_humu395 = ROOT.TGraphErrors(len(humu_indices395), humu_x_plot395, humu_y_plot395, humu_x_plot_err395, humu_y_plot_err395)
gr_humu395.SetMarkerColor(6)
gr_humu395.SetMarkerStyle(21)
gr_humu395.SetMarkerSize(1)

gr_humu789 = ROOT.TGraphErrors(len(humu_indices789), humu_x_plot789, humu_y_plot789, humu_x_plot_err789, humu_y_plot_err789)
gr_humu789.SetMarkerColor(4)
gr_humu789.SetMarkerStyle(21)
gr_humu789.SetMarkerSize(1)

gr_humu1576 = ROOT.TGraphErrors(len(humu_indices1576), humu_x_plot1576, humu_y_plot1576, humu_x_plot_err1576, humu_y_plot_err1576)
gr_humu1576.SetMarkerColor(8)
gr_humu1576.SetMarkerStyle(21)
gr_humu1576.SetMarkerSize(1)

gr_nene = ROOT.TGraphErrors(len(nene_indices), nene_x_plot, nene_y_plot, nene_x_plot_err, nene_y_plot_err)
gr_nene.SetMarkerColor(1)
gr_nene.SetMarkerStyle(22)
gr_nene.SetMarkerSize(1)
fit_nene = gr_nene.Fit("f_lin", "S")
#gr_nene.GetFunction("f_lin").SetParLimits(0, 0, 4000)

gr_nene395 = ROOT.TGraphErrors(len(nene_indices395), nene_x_plot395, nene_y_plot395, nene_x_plot_err395, nene_y_plot_err395)
gr_nene395.SetMarkerColor(6)
gr_nene395.SetMarkerStyle(22)
gr_nene395.SetMarkerSize(1)

gr_nene789 = ROOT.TGraphErrors(len(nene_indices789), nene_x_plot789, nene_y_plot789, nene_x_plot_err789, nene_y_plot_err789)
gr_nene789.SetMarkerColor(4)
gr_nene789.SetMarkerStyle(22)
gr_nene789.SetMarkerSize(1)

gr_nene1576 = ROOT.TGraphErrors(len(nene_indices1576), nene_x_plot1576, nene_y_plot1576, nene_x_plot_err1576, nene_y_plot_err1576)
gr_nene1576.SetMarkerColor(8)
gr_nene1576.SetMarkerStyle(22)
gr_nene1576.SetMarkerSize(1)

mg_FWD = ROOT.TMultiGraph()
mg_FWD.Add(gr_humu)
mg_FWD.Add(gr_humu395)
mg_FWD.Add(gr_humu789)
mg_FWD.Add(gr_humu1576)
mg_FWD.Add(gr_nene)
mg_FWD.Add(gr_nene395)
mg_FWD.Add(gr_nene789)
mg_FWD.Add(gr_nene1576)
#mg_FWD.GetXaxis().SetLimits(0,100000)
mg_FWD.Draw("AP")

l_FWD = ROOT.TLegend(0.1,0.7,0.28,0.6)
l_FWD.AddEntry(gr_humu, "Humu", "Ep")
l_FWD.AddEntry(gr_nene, "Nene", "Ep")
l_FWD.Draw()

###Make histogram of rates using B and T coefficients

B_HER_elepaio = 3.4e2
B_LER_elepaio = 8.7e1
T_HER_elepaio = 0
T_LER_elepaio = 4.4e-5

B_HER_palila =	2.1e3
B_LER_palila =	3.4e2
T_HER_palila =	8.8e-2
T_LER_palila =	3.6e-3

B_HER_tako =	7e1
B_LER_tako =	9.7e1
T_HER_tako =	2.3e-3
T_LER_tako =	0

B_HER_humu =	5.6e4
B_LER_humu =	6.8e2
T_HER_humu =	3.11e-2
T_LER_humu =	1.5e-1

B_HER_nene =	4.1e3
B_LER_nene =	1.3e3
T_HER_nene =	0
T_LER_nene =	1.4e-1



HER_beam_gas_elepaio = B_HER_elepaio * 320 * 1.6e-8 * Z**2
LER_beam_gas_elepaio = B_LER_elepaio * 350 * 6.4e-8 * Z**2
HER_touschek_elepaio = T_HER_elepaio * 320**2/(60*789)
LER_touschek_elepaio = T_LER_elepaio * 350**2/(130*789)

print("Total rates for elepaio are, %s, %s, %s, %s"%(HER_beam_gas_elepaio, LER_beam_gas_elepaio, HER_touschek_elepaio, LER_touschek_elepaio))  

HER_beam_gas_palila = B_HER_palila * 320 * 1.6e-8 * Z**2
LER_beam_gas_palila = B_LER_palila * 350 * 6.4e-8 * Z**2
HER_touschek_palila = T_HER_palila * 320**2/(60*789)
LER_touschek_palila = T_LER_palila * 350**2/(130*789)

print("Total rates for palila are, %s, %s, %s, %s"%(HER_beam_gas_palila, LER_beam_gas_palila, HER_touschek_palila, LER_touschek_palila))

HER_beam_gas_tako = B_HER_tako * 320 * 1.6e-8 * Z**2
LER_beam_gas_tako = B_LER_tako * 350 * 6.4e-8 * Z**2
HER_touschek_tako = T_HER_tako * 320**2/(60*789)
LER_touschek_tako = T_LER_tako * 350**2/(130*789)

print("Total rates for tako are, %s, %s, %s, %s"%(HER_beam_gas_tako, LER_beam_gas_tako, HER_touschek_tako, LER_touschek_tako))

HER_beam_gas_humu = B_HER_humu * 320 * 1.6e-8 * Z**2
LER_beam_gas_humu = B_LER_humu * 350 * 6.4e-8 * Z**2
HER_touschek_humu = T_HER_humu * 320**2/(60*789)
LER_touschek_humu = T_LER_humu * 350**2/(130*789)

print("Total rates for humu are, %s, %s, %s, %s"%(HER_beam_gas_humu, LER_beam_gas_humu, HER_touschek_humu, LER_touschek_humu))

HER_beam_gas_nene = B_HER_nene * 320 * 1.6e-8 * Z**2
LER_beam_gas_nene = B_LER_nene * 350 * 6.4e-8 * Z**2
HER_touschek_nene = T_HER_nene * 320**2/(60*789)
LER_touschek_nene = T_LER_nene * 350**2/(130*789)

print("Total rates for nene are, %s, %s, %s, %s"%(HER_beam_gas_nene, LER_beam_gas_nene, HER_touschek_nene, LER_touschek_nene))
