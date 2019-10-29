### analysis package contains the make_plot(df, Nb, knob) function which generates heuristic plots. IT also includes make_heuristic_dataframe(day,ring,module_id) which is used to generate the dataframe that's passed into the analysis plots

from analysis import *
from os import sys

day = sys.argv[1]
ring = sys.argv[2]
module_id = sys.argv[3]

ROOT.gStyle.SetOptFit(1)

df = make_heuristic_dataframe(day,ring,module_id)

gr = make_plot(df, 789, 0)
gr.SetMarkerStyle(20)
gr.SetMarkerColor(4)

gr1 = make_plot(df, 789, 1)
gr1.SetMarkerStyle(20)
gr1.SetMarkerColor(2)

#gr2 = make_plot(df, 789, 2)
#gr2.SetMarkerStyle(20)
#gr2.SetMarkerColor(1)

gr3 = make_plot(df, 789, -1)
gr3.SetMarkerStyle(20)
gr3.SetMarkerColor(2)

#gr4 = make_plot(df, 789, -2)
#gr4.SetMarkerStyle(20)
#gr4.SetMarkerColor(1)

gr5 = make_plot(df, 1576, 0)
gr5.SetMarkerStyle(22)
gr5.SetMarkerColor(64)

gr6 = make_plot(df, 1576, 1)
gr6.SetMarkerStyle(22)
gr6.SetMarkerColor(96)

#gr7 = make_plot(df, 1576, 2)
#gr7.SetMarkerStyle(22)
#gr7.SetMarkerColor(1)

gr8 = make_plot(df, 1576, -1)
gr8.SetMarkerStyle(22)
gr8.SetMarkerColor(96)

#gr9 = make_plot(df, 1576, -2)
#gr9.SetMarkerStyle(22)
#gr9.SetMarkerColor(1)

mg = ROOT.TMultiGraph()
mg.Add(gr)
mg.Add(gr1)
#mg.Add(gr2)
mg.Add(gr3)
#mg.Add(gr4)
mg.Add(gr5)
mg.Add(gr6)
#mg.Add(gr7)
mg.Add(gr8)
#mg.Add(gr9)
mg.Draw("AP")
mg.Fit("pol1", "FQ")

l = ROOT.TLegend(0.1,0.7,0.28,0.6)
l.AddEntry(gr, "Knob 0, Nb = 789", "Ep")
l.AddEntry(gr1, "Knob +/- 1, Nb = 789", "Ep")
#l.AddEntry(gr2, "Knob +/- 2, Nb = 789", "Ep")
l.AddEntry(gr5, "Knob 0, Nb = 1576", "Ep")
l.AddEntry(gr6, "Knob +/- 1, Nb = 1576", "Ep")
#l.AddEntry(gr7, "Knob +/- 2, Nb = 1576", "Ep")
l.Draw()


