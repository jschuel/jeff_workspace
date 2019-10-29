from root_pandas import read_root
import numpy as np
import pandas as pd
import ROOT
from os import sys

module_id = sys.argv[1]

ROOT.gStyle.SetStatX(0.9);
ROOT.gStyle.SetStatY(0.4);

ch = ROOT.TChain("tracks")

ch.Add("/Users/vahsengrouplaptop/workspace/jeff_workspace/data/6-11_HER/tpc_tools_processed_ntuples/separated_events/whole_study_separated_%s.root"%(module_id))

neutron = ROOT.TCut("recoil_energy < (0.25*length-75) && recoil_energy > (0.0411764*length-64.688)  && recoil_energy > 100")

ch.Draw("length:recoil_energy", neutron, "goff")

energy = ch.GetV2()
length = ch.GetV1()

profile = ROOT.TProfile("profile", "LvE profile %s"%(module_id), 41, 0, 8000, 0 , 40000)
for i in range(0, ch.GetSelectedRows()):
    profile.Fill(energy[i],length[i])

profile.SetMarkerStyle(20)
profile.SetMarkerSize(0.8)
profile.SetLineWidth(2)
profile.SetXTitle("Energy (keV)")
profile.SetYTitle("Length (um)")
profile.SetMaximum(40000)
profile.SetMinimum(0)
profile.Draw()

