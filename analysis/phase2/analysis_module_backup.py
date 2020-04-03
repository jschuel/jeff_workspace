'''Module generates variables for performing single beam background studies'''

import root_pandas as rp
import pandas as pd
import ROOT
import math
import numpy as np
import matplotlib.pyplot as plt
import array

def make_heuristic_plots(df, module_ids):
    ROOT.gStyle.SetOptFit()
    ROOT.gStyle.SetTitleFontSize(0.1)
    f_lin = ROOT.TF1("f_lin", "[0]+[1]*x" )
    l = {}
    tpc_dfs = {}
    gr = {}
    mg = {}
    fit = {}
    x = {}
    y = {}
    x_err = {}
    y_err = {}
    Nbs = [789, 1576]
    for module in module_ids:
        tpc_dfs[module]=pd.DataFrame()
        tpc_dfs[module]['heuristic']= df[module+'_heuristic']
        tpc_dfs[module]['heuristic_err'] =  df[module+'_heuristic_err']
        tpc_dfs[module]['x'] =  df['x']
        tpc_dfs[module]['x_err'] =  df['x_err']
        tpc_dfs[module]['Nb'] = df['Nb']
        tpc_dfs[module] = tpc_dfs[module].loc[tpc_dfs[module]['heuristic']!=0]
        tpc_dfs[module].index=[i for i in range(0,len(tpc_dfs[module]))]
        x[module] = array.array('d', tpc_dfs[module]['x'].to_numpy())
        y[module] = array.array('d', tpc_dfs[module]['heuristic'].to_numpy())
        x_err[module] = array.array('d', tpc_dfs[module]['x_err'].to_numpy())
        y_err[module] = array.array('d', tpc_dfs[module]['heuristic_err'].to_numpy())
        l[module] = ROOT.TLegend()
        mg[module] = ROOT.TMultiGraph()
        gr[module] =  ROOT.TGraphErrors(len(x[module]), x[module], y[module], x_err[module], y_err[module])
        f_lin.SetParLimits(1,0.,1)
        fit[module] = gr[module].Fit("f_lin", "S")
        gr[module].SetTitle('%s'%(module))
        gr[module].SetMarkerStyle(20)
        mg[module].Add(gr[module])
        for Nb in Nbs:
            tpc_dfs[module+'_%s'%(Nb)]=pd.DataFrame()
            tpc_dfs[module+'_%s'%(Nb)] = tpc_dfs[module].loc[tpc_dfs[module]['Nb']== Nb]
            tpc_dfs[module+'_%s'%(Nb)].index = [i for i in range(0,len(tpc_dfs[module+'_%s'%(Nb)]))]
            x[module+'_%s'%(Nb)] = array.array('d', tpc_dfs[module+'_%s'%(Nb)]['x'].to_numpy())
            y[module+'_%s'%(Nb)] = array.array('d', tpc_dfs[module+'_%s'%(Nb)]['heuristic'].to_numpy())
            x_err[module+'_%s'%(Nb)] = array.array('d', tpc_dfs[module+'_%s'%(Nb)]['x_err'].to_numpy())
            y_err[module+'_%s'%(Nb)] = array.array('d', tpc_dfs[module+'_%s'%(Nb)]['heuristic_err'].to_numpy())
            if len(x[module+'_%s'%(Nb)]) > 0:
                gr[module+'_%s'%(Nb)] =  ROOT.TGraphErrors(len(x[module+'_%s'%(Nb)]), x[module+'_%s'%(Nb)], y[module+'_%s'%(Nb)], x_err[module+'_%s'%(Nb)], y_err[module+'_%s'%(Nb)])
                gr[module+'_%s'%(Nb)].SetMarkerStyle(20)
                gr[module+'_%s'%(Nb)].SetName("Nb = %s"%(Nb))
                if Nb == 789:
                    gr[module+'_%s'%(Nb)].SetMarkerColor(4)
                if Nb == 1565:
                    gr[module+'_%s'%(Nb)].SetMarkerColor(8)
                mg[module].Add(gr[module+'_%s'%(Nb)])
                l[module].AddEntry(gr[module+'_%s'%(Nb)],"Nb = %s"%(Nb))
        
        mg[module].SetTitle('%s'%(module))
        mg[module].GetYaxis().SetLabelSize(0.09)
        mg[module].GetXaxis().SetLabelSize(0.09)
        mg[module].SetMinimum(0)
        mg[module].GetXaxis().SetLimits(0,150e3)
    return mg, l, fit

def extract_variables_LUMI(df, module_ids, bin_width):
    bin_dfs = split(df,bin_width) #splits dataframe into chunks of len(bin_width)
    LUMI_df = pd.DataFrame()
    ts = [Bin['ts'].mean() for Bin in bin_dfs]
    dts = [Bin['ts'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    lumi = [Bin['ECL_luminosity'].mean() for Bin in bin_dfs]
    dlumi = [Bin['ECL_luminosity'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    decay = [Bin['decay'].mean() for Bin in bin_dfs]
    #continuous_inj = [Bin['continuous_inj'].mean() for Bin in bin_dfs]
    I_HER = [Bin['SKB_HER_current'].mean() for Bin in bin_dfs]
    dI_HER = [Bin['SKB_HER_current'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    P_HER = [Bin['SKB_HER_pres_avg'].mean() for Bin in bin_dfs]
    dP_HER = [Bin['SKB_HER_pres_avg'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    sy_HER = [Bin['SKB_HER_XRM_sigmay'].mean() for Bin in bin_dfs]
    dsy_HER = [Bin['SKB_HER_XRM_sigmay'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    Nb_HER = [Bin['SKB_HER_NOB'].mean() for Bin in bin_dfs]
    dNb_HER = [Bin['SKB_HER_NOB'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    I_LER = [Bin['SKB_LER_current'].mean() for Bin in bin_dfs]
    dI_LER = [Bin['SKB_LER_current'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    P_LER = [Bin['SKB_LER_pres_avg'].mean() for Bin in bin_dfs]
    dP_LER = [Bin['SKB_LER_pres_avg'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    sy_LER = [Bin['SKB_LER_XRM_sigmay'].mean() for Bin in bin_dfs]
    dsy_LER = [Bin['SKB_LER_XRM_sigmay'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    Nb_LER = [Bin['SKB_LER_NOB'].mean() for Bin in bin_dfs]
    dNb_LER = [Bin['SKB_LER_NOB'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    tpc = {}
    LUMI_df['ts'] = ts
    LUMI_df['dts'] = dts
    LUMI_df['I_HER'] = I_HER
    LUMI_df['dI_HER'] = dI_HER
    LUMI_df['P_HER'] = P_HER
    LUMI_df['dP_HER'] = dP_HER
    LUMI_df['sy_HER'] = sy_HER
    LUMI_df['dsy_HER'] = dsy_HER
    LUMI_df['Nb_HER'] = Nb_HER
    LUMI_df['dNb_HER'] = dNb_HER
    LUMI_df['I_LER'] = I_LER
    LUMI_df['dI_LER'] = dI_LER
    LUMI_df['P_LER'] = P_LER
    LUMI_df['dP_LER'] = dP_LER
    LUMI_df['sy_LER'] = sy_LER
    LUMI_df['dsy_LER'] = dsy_LER
    LUMI_df['Nb_LER'] = Nb_LER
    LUMI_df['dNb_LER'] = dNb_LER
    LUMI_df['lumi'] = lumi
    LUMI_df['dlumi'] = dlumi
    LUMI_df['decay'] = decay
    #LUMI_df['continuous_inj'] = continuous_inj
    for module in module_ids:
        tpc[module+'_mean']=[Bin['%s_neutrons'%(module)].mean() for Bin in bin_dfs]
        tpc[module+'_err']=[Bin['%s_neutrons'%(module)].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
        LUMI_df[module+'_mean'] = tpc[module+'_mean']
        LUMI_df[module+'_err'] = tpc[module+'_err']
    return LUMI_df

def extract_variables(df, module_ids, ring, bin_width): #module_ids is a list of TPCs to include in study
    df = df.loc[df['%s_storage'%(ring)]==1]
    df.index = [i for i in range(0,len(df))]
    bin_dfs = split(df,bin_width) #splits dataframe into chunks of len(bin_width)
    heuristic_df = pd.DataFrame()
    ts = [Bin['ts'].mean() for Bin in bin_dfs]
    dts = [Bin['ts'].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    I = [Bin['SKB_%s_current'%(ring)].mean() for Bin in bin_dfs]
    dI = [Bin['SKB_%s_current'%(ring)].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    P = [Bin['SKB_%s_pres_avg'%(ring)].mean() for Bin in bin_dfs]
    dP = [Bin['SKB_%s_pres_avg'%(ring)].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    sy = [Bin['SKB_%s_XRM_sigmay'%(ring)].mean() for Bin in bin_dfs]
    dsy = [Bin['SKB_%s_XRM_sigmay'%(ring)].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    Nb = [Bin['SKB_%s_NOB'%(ring)].mean() for Bin in bin_dfs]
    dNb = [Bin['SKB_%s_NOB'%(ring)].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
    tpc = {}
    heuristic_df['ts'] = ts
    heuristic_df['dts'] = dts
    heuristic_df['I'] = I
    heuristic_df['dI'] = dI
    heuristic_df['P'] = P
    heuristic_df['dP'] = dP
    heuristic_df['sy'] = sy
    heuristic_df['dsy'] = dsy
    heuristic_df['Nb'] = Nb
    heuristic_df['dNb'] = dNb
    for module in module_ids:
        tpc[module+'_mean']=[Bin['%s_neutrons'%(module)].mean() for Bin in bin_dfs]
        tpc[module+'_err']=[Bin['%s_neutrons'%(module)].std()/math.sqrt(len(Bin)) for Bin in bin_dfs]
        heuristic_df[module+'_mean'] = tpc[module+'_mean']
        heuristic_df[module+'_err'] = tpc[module+'_err']
        heuristic_df[module+'_heuristic'] = heuristic_df[module+'_mean']/(heuristic_df['I']*heuristic_df['P']*2.3**2) #2.3 is Z_eff
        heuristic_df[module+'_heuristic_err']=[0 for i in range(0,len(heuristic_df))] #make errors 0 and repopulate
        for i in range(0,len(heuristic_df)):
            if heuristic_df[module+'_heuristic'][i] == 0:
                heuristic_df[module+'_heuristic_err'][i] = 0
            else:
                heuristic_df[module+'_heuristic_err'][i] = heuristic_df[module+'_heuristic'][i]*math.sqrt((heuristic_df[module+'_err'][i]/heuristic_df[module+'_mean'][i])**2 + (heuristic_df['dI'][i]/heuristic_df['I'][i])**2+(heuristic_df['dP'][i]/heuristic_df['P'][i])**2)
    heuristic_df['x']=heuristic_df['I']/(heuristic_df['P']*heuristic_df['sy']*heuristic_df['Nb']*2.3**2)
    heuristic_df['x_err']=[0 for i in range(0,len(heuristic_df))] #make errors 0 and repopulate
    for i in range(0,len(heuristic_df)):
            if heuristic_df['x'][i] == 0:
                heuristic_df['x_err'][i] = 0
            else:
                heuristic_df['x_err'][i]=heuristic_df['x'][i]*math.sqrt((heuristic_df['dI'][i]/heuristic_df['I'][i])**2+(heuristic_df['dP'][i]/heuristic_df['P'][i])**2+(heuristic_df['dsy'][i]/heuristic_df['sy'][i])**2+(heuristic_df['dNb'][i]/heuristic_df['Nb'][i])**2)
    ### Next line filters out areas where timestamp std error is too long (i.e. areas between multiple fills
    heuristic_df = heuristic_df.loc[heuristic_df['dts']<5]
    heuristic_df.index = [i for i in range(0,len(heuristic_df))]
    #heuristic_df = heuristic_df.loc[heuristic_df['palila_heuristic']<100000]
    #heuristic_df = heuristic_df.loc[heuristic_df['tako_heuristic']>1500]
    heuristic_df.index = [i for i in range(0,len(heuristic_df))]
    ###
    return heuristic_df
        
def split(dfm, bin_width): #splits dataframe into chunks of len(bin_width)
    indices = index_marks(dfm.shape[0], bin_width)
    return np.split(dfm, indices)

def index_marks(nrows, bin_width): #Comes up with index range to use np.split to split dataframe
    return range(bin_width, math.ceil(nrows / bin_width) * bin_width, bin_width)
