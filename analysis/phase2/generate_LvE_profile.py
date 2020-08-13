from root_pandas import read_root
import numpy as np
import pandas as pd
import ROOT
from os import sys

ROOT.gStyle.SetLabelSize(.1, "XY")
ROOT.gStyle.SetTitleSize(.1, "XY")
#ROOT.gStyle.SetTitleOffset(0.6, "Y")
#ROOT.gStyle.SetTitleOffset(0.75, "X")

def get_parameters(module_id):
    ch = ROOT.TChain("tracks")
    ch.Add("~/data/phase2/background_studies/6-11_HER/tpc_tools_processed/separated_ntuples/whole_study_separated_%s.root"%(module_id))
    neutron = ROOT.TCut("recoil_energy < (0.25*length-75) && recoil_energy > (0.0411764*length-64.688)  && recoil_energy > 100")
    ch.Draw("length:recoil_energy", neutron, "goff")
    e = ch.GetV2()
    l = ch.GetV1()
    n = ch.GetSelectedRows()
    energy = [e[i] for i in range(0,n)]
    length = [l[i] for i in range(0,n)]
    return energy, length

c1 = ROOT.TCanvas("c1", "c1", 1400, 1200)
c1.Divide(2,4)

module_id = ["iiwi", "honu", "kohola", "nene", "tako", "humu", "palila", "elepaio"]
positions = ["BWD 198", "BWD 270", "BWD 18", "BWD 90", "FWD 90", "FWD 22", "FWD 270", "FWD 202" ]

ROOT.gStyle.SetStatX(0.9)
ROOT.gStyle.SetStatY(0.5)
ROOT.gStyle.SetStatW(0.25)
ROOT.gStyle.SetStatH(0.25)

profile = []
for i in range(0,len(module_id)):
    c1.cd(i+1)
    profile.append(ROOT.TProfile("profile", "%s"%(positions[i]), 41, 0, 8000, 0 , 40000))
    ROOT.gStyle.SetTitleFontSize(0.1)
    E, L = get_parameters(module_id[i])
    for j in range(0,len(E)):
        profile[i].Fill(E[j],L[j])

    profile[i].SetMarkerStyle(20)
    profile[i].SetMarkerSize(0.5)
    profile[i].SetLineWidth(2)
    #profile[i].SetXTitle("E (keV)")
    #profile[i].SetYTitle("Length (um)")
    profile[i].SetMaximum(40000)
    profile[i].SetMinimum(0)
    profile[i].Draw()
    c1.Update()



