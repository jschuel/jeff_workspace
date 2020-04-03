### analysis package contains the make_plot(df, Nb, knob) function which generates heuristic plots. IT also includes make_heuristic_dataframe(day,ring,module_id) which is used to generate the dataframe that's passed into the analysis plots. It also also includes get_sim_rates which provides simulated rates for various studies

from analysis import *
from os import sys
import array

day = sys.argv[1]
ring = sys.argv[2]

ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetLabelSize(.08, "XY")
ROOT.gStyle.SetOptStat(0)

c1 = ROOT.TCanvas("c1", "c1", 1024, 768)
c1.Divide(2,4)

def get_fit_params():
    module_id = ["iiwi", "honu", "kohola", "nene", "tako", "humu", "palila", "elepaio"]
    positions = [198, 270, 18, 90, 90, 22, 270, 202 ]

    mg = []
    fit = {}
    data_rates = {} #extrapolated rates using B and T from data with simulation beam parameters
    data_rates_BG = {}
    data_rates_Touschek = {}
    data_error = {}

    if ring == "HER":
        sim_params = {'I': 287, 'P': 1.3332e-7, 'sigy': 36, 'Nb': 789} #simulated parameters for extrapolations
    else:
        sim_params = {'I': 341, 'P': 1.3332e-7, 'sigy': 38, 'Nb': 789}

    sim_rates_Coulomb = get_sim_rates("Coulomb", ring)
    sim_rates_Touschek = get_sim_rates("Touschek", ring)
    sim_rates_Brems = get_sim_rates("RBB", ring)
    sim_rates_BG = {}
    sim_rates = {}

    dataMC = {}
    dataMC_BG = {}
    dataMC_Touschek = {}

    for key in sim_rates_Coulomb.keys():
        sim_rates_BG[key] = sim_rates_Brems[key] + sim_rates_Coulomb[key]
        sim_rates[key] = sim_rates_Brems[key] + sim_rates_Coulomb[key] + sim_rates_Touschek[key]

    for i in range(0,len(module_id)):
        mg.append(ROOT.TMultiGraph()) #Creates a MultiGraph object in the list
        
        df = make_heuristic_dataframe(day,ring,module_id[i])
        mg[i] = make_multigraph(df)
        mg[i].SetTitle("%s"%(positions[i]))
        if ring == "HER":
            mg[i].GetXaxis().SetLimits(0,140e3)
            mg[i].SetMinimum(0)
            mg[i].SetMaximum(9000)
        else:
            mg[i].GetXaxis().SetLimits(0,50e3)
            mg[i].SetMinimum(0)
            mg[i].SetMaximum(7000)
        
        f = ROOT.TF1("f","pol1");
        f.SetParLimits(1,0, 10);
        fit[module_id[i]] = mg[i].Fit("f", "SFQ") #adds TPC keys for fit results. Par(0) is B and Par(1) is T
        data_rates[module_id[i]] = fit[module_id[i]].Parameter(0)*sim_params['I']*sim_params['P']*7**2 +  fit[module_id[i]].Parameter(1)*sim_params['I']**2/(sim_params['sigy']*sim_params['Nb'])
        data_error[module_id[i]] = fit[module_id[i]].ParError(0)*sim_params['I']*sim_params['P']*7**2 +  fit[module_id[i]].ParError(1)*sim_params['I']**2/(sim_params['sigy']*sim_params['Nb'])
        data_rates_BG[module_id[i]] = fit[module_id[i]].Parameter(0)*sim_params['I']*sim_params['P']*7**2
        data_rates_Touschek[module_id[i]] = fit[module_id[i]].Parameter(1)*sim_params['I']**2/(sim_params['sigy']*sim_params['Nb'])
        dataMC[module_id[i]] = data_rates[module_id[i]]/sim_rates[module_id[i]]
        dataMC_BG[module_id[i]] = data_rates_BG[module_id[i]]/sim_rates_BG[module_id[i]]
        dataMC_Touschek[module_id[i]] = data_rates_Touschek[module_id[i]]/sim_rates_Touschek[module_id[i]]
    keys = ["data_rates", "data_error", "data_rates_BG", "data_rates_Touschek", "sim_rates", "sim_rates_BG", "sim_rates_Touschek", "dataMC", "dataMC_BG", "dataMC_Touschek", "mg"]
    vals = [data_rates, data_error, data_rates_BG, data_rates_Touschek, sim_rates, sim_rates_BG, sim_rates_Touschek, dataMC, dataMC_BG, dataMC_Touschek, mg]
    params = {}
    for i in range(0,len(keys)):
        params[keys[i]] = vals[i]
    return params

summary_params = get_fit_params()

for i in range(0,8):
    c1.cd(i+1)
    summary_params['mg'][i].Draw("AP")

def create_stack():
    stack = ROOT.THStack("","")
    hs = create_histograms()
    for key in hs.keys():
        hs[key].SetLineWidth(3)
    hs['sim_rates'].SetLineStyle(9)
    hs['sim_rates_BG'].SetLineStyle(9)
    hs['sim_rates_Touschek'].SetLineStyle(9)
    hs['sim_rates'].SetLineColor(1)
    hs['data_rates'].SetLineColor(1)
    hs['sim_rates_BG'].SetLineColor(6)
    hs['data_rates_BG'].SetLineColor(6)
    hs['sim_rates_Touschek'].SetLineColor(7)
    hs['data_rates_Touschek'].SetLineColor(7)
    stack.Add(hs['data_rates'])
    stack.Add(hs['sim_rates'])
    stack.Add(hs['data_rates_BG'])
    stack.Add(hs['sim_rates_BG'])
    stack.Add(hs['data_rates_Touschek'])
    stack.Add(hs['sim_rates_Touschek'])
    #stack.Draw("nostack")
    l = ROOT.TLegend(0.1,0.7,0.48,0.9)
    l.AddEntry(hs['data_rates'], "Total Data", "l")
    l.AddEntry(hs['data_rates_BG'], "Beam-Gas Data", "l")
    l.AddEntry(hs['data_rates_Touschek'], "Touschek Data", "l")
    l.AddEntry(hs['sim_rates'], "Total MC", "l")
    l.AddEntry(hs['sim_rates_BG'], "Beam-Gas MC", "l")
    l.AddEntry(hs['sim_rates_Touschek'], "Touschek MC", "l")
    #l.Draw()
    return stack, l
    
def create_histograms():
    tpcs = ["kohola", "nene", "iiwi", "honu", "humu", "tako", "elepaio", "palila"]
    phis = ["BWD 18", "BWD 90", "BWD 198", "BWD 270", "FWD 22", "FWD 90", "FWD 202", "FWD 270"]
    params = get_fit_params()
    hists = {}
    xs = [i for i in range(1,9)]
    for key in params.keys():
        if key != "mg":
            hists[key] = ROOT.TH1F("", "", 8, -0.5, 7.5)
            for i in range(0,8):
                hists[key].SetBinContent(xs[i], params[key][tpcs[i]])
                if key == "data_error":
                    hists["data_rates"].SetBinError(xs[i], params[key][tpcs[i]])
                hists[key].GetXaxis().SetBinLabel(i+1, phis[i])
    return hists


'''
hs1 = ROOT.THStack("hs1","Data MC Comparison Beam-Gas Sensitivities; TPC Number;  B-G Sensitivity [Hz/(mA Pa)]")

h_B_data = ROOT.TH1F("h_B_data", "", 8, -0.5, 7.5)
h_B_data.SetStats(0)
h_B_data.SetLineColor(1)
h_B_data.SetLineWidth(8)
h_B_data.SetMarkerStyle(20)
h_B_data.SetMarkerSize(2)
for i in range(0,8):
    h_B_data.SetBinContent(x[i], bg_sensitivities_data[i])
    h_B_data.SetBinError(x[i], bg_sensitivities_data_err[i])
h_B_data.Draw("P")

h_B_sim = ROOT.TH1F("h_B_sim", "", 8, -0.5, 7.5)
h_B_sim.SetStats(0)
h_B_sim.SetLineColor(4)
h_B_sim.SetLineWidth(8)
h_B_sim.SetLineStyle(10)
for i in range(0,8):
    h_B_sim.SetBinContent(x[i], bg_sensitivities_sim[i])
h_B_sim.Draw("P")

hs1.Add(h_B_data)
hs1.Add(h_B_sim)
hs1.Draw("nostack")

l1 = ROOT.TLegend(0.1,0.7,0.48,0.9)
l1.AddEntry(h_B_data, "B (Data)", "l")
l1.AddEntry(h_B_sim, "B (MC)", "l")
l1.Draw()

c3.cd(2)

hs2 = ROOT.THStack("hs2","Data MC Comparison Touschek Sensitivities; TPC Number;  Touschek Sensitivity [(Hz um)/(mA^2)]")

h_T_data = ROOT.TH1F("h_T_data", "", 8, -0.5, 7.5)
h_T_data.SetStats(0)
h_T_data.SetLineColor(1)
h_T_data.SetLineWidth(8)
h_T_data.SetMarkerStyle(20)
h_T_data.SetMarkerSize(2)
for i in range(0,8):
    h_T_data.SetBinContent(x[i], touschek_sensitivities_data[i])
    h_T_data.SetBinError(x[i], touschek_sensitivities_data_err[i])
h_T_data.Draw("P")

h_T_sim = ROOT.TH1F("h_T_sim", "", 8, -0.5, 7.5)
h_T_sim.SetStats(0)
h_T_sim.SetLineColor(4)
h_T_sim.SetLineWidth(8)
h_T_sim.SetLineStyle(10)
for i in range(0,8):
    h_T_sim.SetBinContent(x[i], touschek_sensitivities_sim[i])
h_T_sim.Draw("P")

hs2.Add(h_T_data)
hs2.Add(h_T_sim)
hs2.Draw("nostack")

l2 = ROOT.TLegend(0.1,0.7,0.48,0.9)
l2.AddEntry(h_T_data, "T (Data)", "l")
l2.AddEntry(h_T_sim, "T (MC)", "l")
l2.Draw()
'''
