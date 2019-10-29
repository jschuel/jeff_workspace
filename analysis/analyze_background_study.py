### analysis package contains the make_plot(df, Nb, knob) function which generates heuristic plots. IT also includes make_heuristic_dataframe(day,ring,module_id) which is used to generate the dataframe that's passed into the analysis plots

from analysis import *
from os import sys

day = sys.argv[1]
ring = sys.argv[2]

ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetLabelSize(.08, "XY")

c1 = ROOT.TCanvas("c1", "c1", 1024, 768)
c1.Divide(2,4)

module_id = ["iiwi", "honu", "kohola", "nene", "tako", "humu", "palila", "elepaio"]

mg = []

for i in range(0,len(module_id)):
    mg.append(ROOT.TMultiGraph()) #Creates a MultiGraph object in the list
    c1.cd(i+1)
    df = make_heuristic_dataframe(day,ring,module_id[i])
    mg[i] = make_multigraph(df)
    if ring == "HER":
        mg[i].GetXaxis().SetLimits(0,140e3)
        mg[i].SetMinimum(0)
        mg[i].SetMaximum(9000)
    else:
        mg[i].GetXaxis().SetLimits(0,50e3)
        mg[i].SetMinimum(0)
        mg[i].SetMaximum(7000)
    mg[i].Draw("AP")
    mg[i].Fit("pol1", "FQ")


