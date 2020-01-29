import ROOT
import numpy as np
import matplotlib as plt
import sys
import math
import pandas as pd
from root_pandas import read_root
from analysis import *

###Parameters for startup

date = sys.argv[1]
ring = sys.argv[2]


### GENERATE DATA  ###


df = read_root("/Users/vahsengrouplaptop/data/phase2/combined_SKB_TPC_ntuples/June_%s_%s_all_pressures.root"%(date, ring))
df.index = [i for i in range(0,len(df))]
index = [i for i in range (3693, len(df))]
df = df.drop(index, axis=0)
index = [i for i in range(0, 151)] + [i for i in range(1250, 1420)] + [i for i in range(2750, 2850)]
df = df.drop(index, axis=0)
df.index = [i for i in range(0,len(df))]
df.columns = ['ts', 'I', 'P_avg', 'sigy', 'P_D1', 'P_D2', 'P_D3', 'P_D4', 'P_D5',
       'P_D6', 'P_D7', 'P_D8', 'P_D9', 'P_D10', 'P_D11', 'P_D12', 'Nb',
       'knob', 'iiwi_neutrons', 'honu_neutrons', 'kohola_neutrons',
       'nene_neutrons', 'tako_neutrons', 'humu_neutrons', 'palila_neutrons',
       'elepaio_neutrons']
'''
#chunk makes all pressures equal to D02 to see the effects of increasing pressure in a sector when the pressure is initially uniform around the ring
for i in range(0,12):
    df['P_D%s'%(i+1)] = df['P_D2']
df['P_avg'] = df['P_D2']
'''
#df['P_D2']= 5*df['P_D2']
#df['P_avg'] = 1.333333333*df['P_avg']

He_rates_C = pd.DataFrame.from_dict(get_pressure_reweight_sim_rates("Coulomb", ring, "He"))
C_rates_C = pd.DataFrame.from_dict(get_pressure_reweight_sim_rates("Coulomb", ring, "C"))
O_rates_C = pd.DataFrame.from_dict(get_pressure_reweight_sim_rates("Coulomb", ring, "O"))
He_rates_B = pd.DataFrame.from_dict(get_pressure_reweight_sim_rates("Brems", ring, "He"))
C_rates_B = pd.DataFrame.from_dict(get_pressure_reweight_sim_rates("Brems", ring, "C"))
O_rates_B = pd.DataFrame.from_dict(get_pressure_reweight_sim_rates("Brems", ring, "O"))
He_rates_T = pd.DataFrame.from_dict(get_pressure_reweight_sim_rates("Touschek", ring, "He"))
C_rates_T = pd.DataFrame.from_dict(get_pressure_reweight_sim_rates("Touschek", ring, "C"))
O_rates_T = pd.DataFrame.from_dict(get_pressure_reweight_sim_rates("Touschek", ring, "O"))

sim_rates_BG = He_rates_C + C_rates_C + O_rates_C + He_rates_B + C_rates_B + O_rates_B
    
#for col in sim_rates_BG.columns:
#    sim_rates_BG[col][1] = 5*sim_rates_BG[col][1]

sim_rates_T = He_rates_T + C_rates_T + O_rates_T


sim_rates = sim_rates_BG + sim_rates_T

###Parameters from simulation###
    
if ring == "HER":
    sim_params = {'I': 287, 'P': 1.3332e-7, 'sigy': 36, 'Nb': 789} #simulated parameters for extrapolations
else:
    sim_params = {'I': 341, 'P': 1.3332e-7, 'sigy': 38, 'Nb': 789}
I = sim_params['I']
P = sim_params['P']
sigy = sim_params['sigy']
Nb = sim_params['Nb']
Ze = 7

###Base pressures (Comment out if you don't want to modify to pressure seen by beam)###

#P_base = {'P_HER_D1': 1.85e-08, 'P_HER_D2': 1.88e-08, 'P_HER_D3': 1e-08, 'P_HER_D4': 1.62e-08, 'P_HER_D5': 1.54e-08, 'P_HER_D6': 1e-08, 'P_HER_D7': 1e-08, 'P_HER_D8': 1.12e-08, 'P_HER_D9': 1.32e-08, 'P_HER_D10': 1e-08, 'P_HER_D11': 1e-08, 'P_HER_D12': 1.1e-08, 'P_LER_D1': 3.5e-08, 'P_LER_D2': 3.33e-08, 'P_LER_D3': 1e-08, 'P_LER_D4': 2e-08, 'P_LER_D5': 2.36e-08, 'P_LER_D6': 1.02e-08, 'P_LER_D7': 4.65e-08, 'P_LER_D8': 2.52e-08, 'P_LER_D9': 1e-08, 'P_LER_D10': 2.42e-08, 'P_LER_D11': 1.92e-08, 'P_LER_D12': 1.04e-08} #from may beam off period

P_base = {'P_HER_D1': 1.104e-08, 'P_HER_D2': 1.0404e-08, 'P_HER_D3': 1e-08, 'P_HER_D4': 1.273e-08, 'P_HER_D5': 1e-08, 'P_HER_D6': 1e-08, 'P_HER_D7': 1e-08, 'P_HER_D8': 1.285e-08, 'P_HER_D9': 1.329e-08, 'P_HER_D10': 1.057e-08, 'P_HER_D11': 1e-08, 'P_HER_D12': 1.096e-08, 'P_LER_D1': 1.863e-08, 'P_LER_D2': 1.177e-08, 'P_LER_D3': 1.022e-08, 'P_LER_D4': 1.256e-08, 'P_LER_D5': 2.408e-08, 'P_LER_D6': 1e-08, 'P_LER_D7': 2.024e-08, 'P_LER_D8': 2.296e-08, 'P_LER_D9': 1e-08, 'P_LER_D10': 1.517e-08, 'P_LER_D11': 1.315e-08, 'P_LER_D12': 1.048e-08} #from post phase 2 10 day beam off period

''' #below is for method 2 of obtaining base pressures (y-int of fit to P vs I)
P_base = []
for i in range(1,13):
    fit = np.polyfit(df['I'], df['P_D%s'%(i)],1)
    P_base.append(fit[1])
'''

###Create weighted pressures
weights = pd.DataFrame()
pseudo_rates = pd.DataFrame()
reweighted_pressures = pd.DataFrame()

for col in sim_rates_BG.columns:
    weights[col] = sim_rates_BG[col]/sum(sim_rates_BG[col])
    pseudo_rates[col] = (1/sim_params['P'])*sum([sim_rates_BG[col][i]*df['P_D%s'%(i+1)] for i in range(0,len(sim_rates_BG[col]))])+sum(sim_rates_T[col])
    reweighted_pressures[col] = sum([weights[col][i]*df['P_D%s'%(i+1)] for i in range(0,len(sim_rates_BG[col]))])
''' #Chunk below keeps pseudo rates the same while modifying pressure to "beam" pressure
for col in sim_rates_BG.columns:
    weights[col] = sim_rates_BG[col]/sum(sim_rates_BG[col])
    pseudo_rates[col] = (1/sim_params['P'])*sum([sim_rates_BG[col][i]*df['P_D%s'%(i+1)] for i in range(0,len(sim_rates_BG[col]))])+sum(sim_rates_T[col])
for i in range(1,13):
    #df['P_D%s'%i] = 3*df['P_D%s'%i] - 2*P_base[i-1] #method 2
    if ring == "HER": #method 1
        df['P_D%s'%i] = 3*df['P_D%s'%i] - 2*P_base['P_HER_D%s'%i]
    else:
        df['P_D%s'%i] = 3*df['P_D%s'%i] - 2*P_base['P_LER_D%s'%i]
for col in sim_rates_BG.columns:
    reweighted_pressures[col] = sum([weights[col][i]*df['P_D%s'%(i+1)] for i in range(0,len(sim_rates_BG[col]))])
'''
###Define ideal heuristic parameters

x_ideal = np.array(([0, I/(P*sigy*Nb*Ze**2)]))
y_ideal = pd.DataFrame()

for col in sim_rates_BG.columns:
    y_ideal[col] = [np.sum(sim_rates_BG[col].to_numpy())/(I*P*Ze**2),np.sum(sim_rates[col].to_numpy())/(I*P*Ze**2)]



y_obs = pd.DataFrame()
x_obs = pd.DataFrame()
x_obs_reweight = pd.DataFrame()
y_obs_reweight = pd.DataFrame()


for col in sim_rates_BG.columns:
    y_obs[col] = [pseudo_rates[col][j]/(I*df['P_avg'][j]*Ze**2) for j in range(0,len(df))]
    x_obs[col] = [I/(df['P_avg'][j]*sigy*Nb*Ze**2) for j in range(0,len(df))]
    x_obs_reweight[col] = [I/(reweighted_pressures[col][j]*sigy*Nb*Ze**2) for j in range(0,len(df))]
    y_obs_reweight[col] = [(pseudo_rates[col][j])/(I*reweighted_pressures[col][j]*Ze**2) for j in range(0,len(df))]

c1 = ROOT.TCanvas("c1", "c1", 1024, 768)
c1.Divide(2,4)

module_id = [col for col in sim_rates_BG.columns]

fit_avg = {} #dictionary of fit parameters for ring averaged pressure study
fit_reweight = {} #dictionary of fit parameters for reweighted study
fit_ideal = {}

mg = [ROOT.TMultiGraph() for i in range(0,len(module_id))]

for i in range(0,len(module_id)):
    c1.cd(i+1)
    gr_ideal = ROOT.TGraph(2, x_ideal, y_ideal[module_id[i]].to_numpy())
    gr_ideal.SetMarkerStyle(20)
    gr_ideal.SetMarkerColor(4)
    gr_ideal.SetMarkerSize(0)
    gr_ideal.SetLineColor(4)
    gr_ideal.GetXaxis().SetLimits(0,500000)
    f_ideal = ROOT.TF1("f_ideal", "[0]+[1]*x" )
    r_ideal = gr_ideal.Fit("f_ideal", "S")
    fit_ideal[module_id[i]] = gr_ideal.Fit("f_ideal", "S")
    gr_ideal.GetFunction("f_ideal").SetLineColor(4)


    gr_obs = ROOT.TGraph(len(x_obs), x_obs[module_id[i]].to_numpy(), y_obs[module_id[i]].to_numpy())
    gr_obs.SetMarkerStyle(20)
    gr_obs.SetMarkerColor(4)
    gr_obs.SetMarkerSize(0.5)
    gr_obs.GetXaxis().SetLimits(0,500000)
    f_obs = ROOT.TF1("f_obs", "[0]+[1]*x" )
    fit_avg[module_id[i]] = gr_obs.Fit("f_obs", "S")
    gr_obs.GetFunction("f_obs").SetLineColor(2)

    gr_reweight = ROOT.TGraph(len(x_obs_reweight), x_obs_reweight[module_id[i]].to_numpy(), y_obs_reweight[module_id[i]].to_numpy())
    gr_reweight.SetMarkerStyle(20)
    gr_reweight.SetMarkerColor(1)
    gr_reweight.SetMarkerSize(0.5)
    f_reweight = ROOT.TF1("f_reweight", "[0]+[1]*x" )
    fit_reweight[module_id[i]] = gr_reweight.Fit("f_reweight", "S")

    mg[i] = ROOT.TMultiGraph()
    mg[i].Add(gr_ideal)
    #mg[i].Add(gr_obs)
    mg[i].Add(gr_reweight)
    mg[i].GetXaxis().SetLabelSize(0.08)
    mg[i].GetYaxis().SetLabelSize(0.08)
    c1.Update()
    mg[i].SetMinimum(0)
    #mg[i].SetMaximum(6000)
    mg[i].GetXaxis().SetLimits(0,15000)
    mg[i].Draw("AP")
    
'''    
l = ROOT.TLegend(0.1,0.7,0.28,0.6)
l.AddEntry(gr_ideal, "Predicted Combined heuristic at simulated beam conditions", "l")
l.AddEntry(gr_obs, "Pseudo-rates scaled to spatially averaged pressure with simulated beam currents", "p")
l.AddEntry(gr_reweight, "Pseudo-rates scaled to reweighted pressure with simulated beam currents", "p")
l.Draw()
'''

'''
bg_ratio = {}
T_ratio = {}
for i in range(0,len(module_id)):
    bg_ratio[module_id[i]] = '%.3g' % (fit_reweight[module_id[i]].Parameter(0)/(fit_ideal[module_id[i]].Parameter(0)))
    T_ratio[module_id[i]] = '%.3g' % (fit_reweight[module_id[i]].Parameter(1)/(fit_ideal[module_id[i]].Parameter(1)))
'''

bg_params = {}
T_params = {}
for i in range(0,len(module_id)):
    bg_params[module_id[i]] = '%.3g' % (fit_reweight[module_id[i]].Parameter(0))
    T_params[module_id[i]] = '%.2g' % (fit_reweight[module_id[i]].Parameter(1))
