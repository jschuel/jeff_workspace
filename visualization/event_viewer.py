import ROOT
import numpy as np
import array
import root_pandas as rp
from os import sys

def view_event(entry, ch): 
    c1 = ROOT.TCanvas("c1", "c1", 1600, 1000)
    c1.Divide(2,3)
    entry_cut = ROOT.TCut("event_number == %s"%(entry))
    h1 = ROOT.TH2F("h1","2D event display", 80, 0, 80, 336, 0, 336)
    h2 = ROOT.TH1F("h2","BCID distribution", 100, 0, 100)
    h3 = ROOT.TH1F("h3","TOT distribution", 14,0,14)
    ch.SetMarkerStyle(20)

    c1.cd(1)
    ch.Draw("row:column:tot>>h1", entry_cut ,"COLZ")

    c1.cd(2)
    ch.Draw("BCID:row:column:tot", entry_cut)

    c1.cd(3)
    ch.Draw("BCID>>h2", entry_cut)

    c1.cd(4)
    ch.Draw("tot>>h3", entry_cut)

    c1.cd(5)
    ch.Draw("track_energy:length", entry_cut, "goff");

    single_length = ch.GetV2()
    single_energy = ch.GetV1()

    gr = ROOT.TGraph(1, single_length, single_energy) #to be included in multigraph: plots event in red
    gr.SetMarkerStyle(20)
    gr.SetMarkerColor(2)

    c1.cd(6)
    ch.Draw("track_energy:length")
    length = ch.GetV2()
    energy = ch.GetV1()
    n = ch.GetSelectedRows()
  
    gr1 = ROOT.TGraph(n, length, energy)
    gr1.SetMarkerStyle(20)
    gr1.SetMarkerSize(0.1)
    gr1.SetMarkerColor(4)

    mg = ROOT.TMultiGraph()
    mg.Add(gr1)
    mg.Add(gr)
    mg.SetTitle("Energy vs. Length; Length(um); Recoil Energy(keV)")
    #mg.GetXaxis().SetLimits(0,3000)
    #mg.SetMaximum(800)
    mg.SetMinimum(-2)
    mg.Draw("AP")
    input("Press enter to continue…")
    

def get_TChain(f):
    ch = ROOT.TChain("data")
    ch.Add(f)
    return ch

module_id = sys.argv[1]
run_num = sys.argv[2]

#f = "~/%s_%s.root"%(run_num, module_id)
f = "~/data/phase3/spring_2020/05-09-20/tpc_root_files/%s_%s.root"%(run_num, module_id)
alpha = ROOT.TCut("hitside_col_min == 1 && hitside_col_max == 1 && hitside_row_min == 0 && hitside_row_max == 0 && track_energy > 50")
xray = ROOT.TCut("hitside_col_min == 0 && hitside_col_max == 0 && hitside_row_min == 0 && hitside_row_max == 0 && track_energy < 10")
neutron = ROOT.TCut("track_energy < (0.7*length-75) && track_energy > (0.015*length-65)  && track_energy > 20 && hitside_col_min == 0 && hitside_col_max == 0 && hitside_row_min == 0 && hitside_row_max == 0")
ch = get_TChain(f)
