from root_pandas import read_root
import numpy as np
import pandas as pd
import ROOT
import array
from os import sys

ROOT.gStyle.SetLabelSize(.06, "XY")
ROOT.gStyle.SetTitleSize(.06, "XY")
ROOT.gStyle.SetTitleOffset(0.75, "Y")
ROOT.gStyle.SetTitleOffset(0.85, "X")

def get_parameters(module_id):
    ch = ROOT.TChain("data")
    neutrons = ROOT.TChain("data")
    ch.Add("/Users/vahsengrouplaptop/data/phase2/background_studies/7-12_lumi/%s_all.root"%(module_id))
    neutrons.Add("/Users/vahsengrouplaptop/data/phase2/background_studies/7-12_lumi/%s_all.root"%(module_id))
    neutron = ROOT.TCut("track_energy < (0.5*length-75) && track_energy > (0.04*length-65)  && track_energy > 98.375 && hitside_col_min == 0 && hitside_col_max == 0 && hitside_row_min == 0 && hitside_row_max == 0")
    neutrons.Draw("length:track_energy", neutron, "goff")                                                                                                                                        
    ch.Draw("length:track_energy", "hitside_col_min == 0 && hitside_row_min == 0 && hitside_col_max == 0 && hitside_row_max == 0", "goff")
    #ch.Add("/Users/vahsengrouplaptop/data/phase3/phase3_background_root/tpc_tools/Dec_7_%s_phase3.root"%(module_id))
    #neutrons.Add("/Users/vahsengrouplaptop/data/phase3/phase3_background_root/tpc_tools/Dec_7_%s_phase3.root"%(module_id))
    #neutron = ROOT.TCut("(recoil_energy/1000) < (0.5*length-75) && (recoil_energy/1000) > (0.0411764*length-64.688)  && (recoil_energy/1000) > 100 && hitside_top == 0 && hitside_bottom == 0 && hitside_source == 0 && hitside_antisource == 0 && num_clusters == 1")
    #neutrons.Draw("length:(recoil_energy/1000)", neutron, "goff")
    #ch.Draw("length:(recoil_energy/1000)", "hitside_top == 0 && hitside_bottom == 0 && hitside_source == 0 && hitside_antisource == 0", "goff")
    e = ch.GetV2()
    l = ch.GetV1()
    n = ch.GetSelectedRows()
    e_n = neutrons.GetV2()
    l_n = neutrons.GetV1()
    n_n = neutrons.GetSelectedRows()
    energy = array.array('d', [e[i] for i in range(0,n)])
    length = array.array('d', [l[i] for i in range(0,n)])
    neutron_energy = array.array('d', [e_n[i] for i in range(0,n_n)])
    neutron_length = array.array('d', [l_n[i] for i in range(0,n_n)])
    return energy, length, neutron_energy, neutron_length

#module_id_bwd = ["tako", "palila", "elepaio"]
#module_id_fwd = ["iiwi", "humu", "nene"]
module_id_bwd = ["iiwi", "honu", "kohola", "nene"]
module_id_fwd = ["tako", "humu", "palila", "elepaio"]
gr_bwd = []
gr_neutron_bwd = []
gr_fwd = []
gr_neutron_fwd = []
mg = ROOT.TMultiGraph()
for i in range(0,len(module_id_bwd)):
    energy_bwd, length_bwd, neutron_energy_bwd, neutron_length_bwd = get_parameters(module_id_bwd[i])
    gr_bwd.append(ROOT.TGraph(len(energy_bwd), length_bwd, energy_bwd))
    gr_bwd[i].SetMarkerStyle(20)
    gr_bwd[i].SetMarkerSize(0.1)
    gr_bwd[i].SetMarkerColor(1)
    gr_neutron_bwd.append(ROOT.TGraph(len(neutron_energy_bwd), neutron_length_bwd, neutron_energy_bwd))
    gr_neutron_bwd[i].SetMarkerStyle(20)
    gr_neutron_bwd[i].SetMarkerSize(0.1)
    gr_neutron_bwd[i].SetMarkerColor(4)
    mg.Add(gr_bwd[i])
for i in range(0,len(module_id_fwd)):
    energy_fwd, length_fwd, neutron_energy_fwd, neutron_length_fwd = get_parameters(module_id_fwd[i])
    gr_fwd.append(ROOT.TGraph(len(energy_fwd), length_fwd, energy_fwd))
    gr_fwd[i].SetMarkerStyle(20)
    gr_fwd[i].SetMarkerSize(0.1)
    gr_fwd[i].SetMarkerColor(15)
    gr_neutron_fwd.append(ROOT.TGraph(len(neutron_energy_fwd), neutron_length_fwd, neutron_energy_fwd))
    gr_neutron_fwd[i].SetMarkerStyle(20)
    gr_neutron_fwd[i].SetMarkerSize(0.1)
    gr_neutron_fwd[i].SetMarkerColor(8)
    mg.Add(gr_fwd[i])
for i in range(0,len(module_id_bwd)):
    mg.Add(gr_neutron_bwd[i])
for i in range(0,len(module_id_fwd)):
    mg.Add(gr_neutron_fwd[i])
mg.SetMaximum(3000)
mg.SetMinimum(-20)
mg.Draw("AP")

