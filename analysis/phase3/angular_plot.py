import ROOT
import pandas as pd
import root_pandas as rp
import numpy as np
import array
import math
from os import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

module = sys.argv[1]

def make_plot(ring, module_id):
    if ring == "LUMI":
        df = rp.read_root("/Users/vahsengrouplaptop/data/phase3/PVM/Dec_8_LUMI_updated.root")
        index = df.loc[(df['ECL_LUM_MON_ACC']>7000) & (df['continuous_inj'] == 1)].index.to_numpy()
        t_low = df['ts'][index].min()
        t_high = df['ts'][index].max()
        ch = ROOT.TChain("tracks")
        ch1 = ROOT.TChain("tracks")
        ch.Add("/Users/vahsengrouplaptop/data/phase3/phase3_background_root/tpc_tools/Dec_8_%s_phase3.root"%(module_id))
        ch1.Add("/Users/vahsengrouplaptop/data/phase3/phase3_background_root/tpc_tools/Dec_8_%s_phase3.root"%(module_id))
        ch.Draw("phi:theta:recoil_energy", "ts >%s && ts < %s && (recoil_energy/1000) < (0.5*length-75) && (recoil_energy/1000) > (0.0411764*length-64.688)  && (recoil_energy/1000) > 100 && hitside_top == 0 && hitside_bottom == 0 && hitside_source == 0 && hitside_antisource == 0 && num_clusters == 1"%(t_low,t_high), "goff")
        theta = ch.GetV1()
        phi = ch.GetV2()
        energy = ch.GetV3()
        ch1.Draw("head_charge:tail_charge", "ts >%s && ts < %s && (recoil_energy/1000) < (0.5*length-75) && (recoil_energy/1000) > (0.0411764*length-64.688)  && (recoil_energy/1000) > 100 && hitside_top == 0 && hitside_bottom == 0 && hitside_source == 0 && hitside_antisource == 0 && num_clusters == 1"%(t_low,t_high), "goff")
        head_charge = ch1.GetV1()
        tail_charge = ch1.GetV2()
        vec = ROOT.TVector3()
        scatter = ROOT.TH2F("scatter", "High Lumi Phi vs. theta; theta; phi",20,-1.1,1.1,20,-185,185)
        theta_dist = ROOT.TH1F("theta_dist", "High Lumi Theta Distribution; theta", 20, -5, 185)
        phi_dist = ROOT.TH1F("phi_dist", "High Lumi Phi Distribution; phi", 20, -185, 185)
        n = ch.GetSelectedRows()
        for i in range(0,n):
            vec.SetMagThetaPhi(1,theta[i],phi[i])
            #if(vec.Theta() > np.pi/2):
            if (head_charge[i] > tail_charge[i]):
                vec = -1*vec
            theta[i] = vec.Theta()*180/np.pi
            phi[i] = vec.Phi()*180/np.pi
            scatter.Fill(np.cos(theta[i]*np.pi/180),phi[i])
            #scatter.Fill(theta[i],phi[i])
            theta_dist.Fill(theta[i])
            theta_dist.SetMinimum(0)
            theta_dist.SetLineWidth(3)
            phi_dist.Fill(phi[i])
            phi_dist.SetMinimum(0)
            phi_dist.SetLineWidth(3)
            
    else:
        df_tpc = rp.read_root("/Users/vahsengrouplaptop/data/phase3/phase3_background_root/tpc_tools/Dec_7_%s_phase3.root"%(module_id),"tracks")
        ch = ROOT.TChain("tracks")
        ch.Add("/Users/vahsengrouplaptop/data/phase3/phase3_background_root/tpc_tools/Dec_7_%s_phase3.root"%(module_id))
        ch.Draw("ts", "(recoil_energy/1000) < (0.5*length-75) && (recoil_energy/1000) > (0.0411764*length-64.688)  && (recoil_energy/1000) > 100 && hitside_top == 0 && hitside_bottom == 0 && hitside_source == 0 && hitside_antisource == 0 && num_clusters == 1", "goff")
        time = ch.GetV1()
        ts = [time[i] for i in range(0,ch.GetSelectedRows())]
        df_tpc = df_tpc[df_tpc['ts'].isin(ts)]
        df_tpc.index = [i for i in range(0,len(df_tpc))]
        df_tpc['ts'] = [int(val) for val in df_tpc['ts']]
        df = rp.read_root("/Users/vahsengrouplaptop/data/phase3/PVM/Dec_7_%s_updated.root"%(ring))
        ts_index = df.loc[df['Storage_Flag'] == 1].index.to_numpy()
        #ts_index = df.loc[df['Storage_Flag'] < 2].index.to_numpy()
        timestamps = [int(ts) for ts in df['ts'][ts_index]]
        df_tpc = df_tpc[df_tpc['ts'].isin(timestamps)]
        df_tpc.index = [i for i in range(0,len(df_tpc))]
        theta = array.array('d', df_tpc['theta'])
        theta_plt = array.array('d', df_tpc['theta'])
        phi = array.array('d', df_tpc['phi'])
        phi_plt = array.array('d', df_tpc['phi'])
        energy = array.array('d', df_tpc['recoil_energy'])
        head_charge = array.array('d', df_tpc['head_charge'])
        tail_charge = array.array('d', df_tpc['tail_charge'])
        scatter = ROOT.TH2F("scatter_%s"%(ring), "%s Beam Storage Phi vs. theta; theta; phi"%(ring),20,-1.1,1.1,20,-185,185)
        theta_dist = ROOT.TH1F("theta_dist_%s"%(ring), "%s Beam Storage Theta Distribution; theta"%(ring), 20, -5, 185)
        phi_dist = ROOT.TH1F("phi_dist_%s"%(ring), "%s Beam Storage Phi Distribution; phi"%(ring), 20, -185, 185)
        n = len(df_tpc)
        vec = ROOT.TVector3()
        for i in range(0,n):
            vec.SetMagThetaPhi(1,theta[i],phi[i])
            #if(vec.Theta() > np.pi/2):
            if (head_charge[i] > tail_charge[i]):
                vec = -1*vec
            theta[i] = vec.Theta()*180/np.pi
            theta_plt[i] = vec.Theta()
            phi[i] = vec.Phi()*180/np.pi
            phi_plt[i] = vec.Phi()
            scatter.Fill(np.cos(theta[i]*np.pi/180),phi[i])
            #scatter.Fill(theta[i],phi[i])
            theta_dist.Fill(theta[i])
            theta_dist.SetMinimum(0)
            theta_dist.SetLineWidth(3)
            phi_dist.Fill(phi[i])
            phi_dist.SetMinimum(0)
            phi_dist.SetLineWidth(3)
        ''' Spherical projection plot
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = np.cos(u)*np.sin(v)
        y = np.sin(u)*np.sin(v)
        z = np.cos(v)
        ax.plot_wireframe(x, y, z, color="r")
        ax.scatter(np.sin(theta_plt)*np.cos(phi_plt), np.sin(theta_plt)*np.sin(phi_plt), np.cos(theta_plt), s = 1)
        plt.show()
        '''
    if ring == "LUMI":
        return scatter, theta_dist, phi_dist
    else:
        return scatter, theta_dist, phi_dist, df_tpc, df

if module != "iiwi" and module != "nene" and module != "humu":
    colz, theta, phi = make_plot("LUMI", module)
    colz_LER, theta_LER, phi_LER = make_plot("LER", module)
    colz_HER, theta_HER, phi_HER = make_plot("HER", module)
    c1 = ROOT.TCanvas('c1', '%s'%(module), 800, 600)
    c1.Divide(3,3)
    c1.cd(1)
    colz.Draw("COLZ")
    c1.cd(2)
    colz_LER.Draw("COLZ")
    c1.cd(3)
    colz_HER.Draw("COLZ")
    c1.cd(4)
    theta.Draw()
    c1.cd(5)
    theta_LER.Draw()
    c1.cd(6)
    theta_HER.Draw()
    c1.cd(7)
    phi.Draw()
    c1.cd(8)
    phi_LER.Draw()
    c1.cd(9)
    phi_HER.Draw()
else:
    colz, theta, phi = make_plot("LUMI", module)
    colz_LER, theta_LER, phi_LER, df_tpc, df_SKB = make_plot("LER", module)
    c1 = ROOT.TCanvas('c1', '%s'%(module), 800, 600)
    c1.Divide(2,3)
    c1.cd(1)
    colz.Draw("COLZ")
    c1.cd(2)
    colz_LER.Draw("COLZ")
    c1.cd(3)
    theta.Draw()
    c1.cd(4)
    theta_LER.Draw()
    c1.cd(5)
    phi.Draw()
    c1.cd(6)
    phi_LER.Draw()
