import root_pandas as rp
import pandas as pd
import ROOT
import math
import numpy as np
import matplotlib.pyplot as plt
import array
from analysis_module import extract_variables
from analysis_module import extract_variables_LUMI
from analysis_module import make_heuristic_plots
from os import sys

month = sys.argv[1]
day = sys.argv[2]

input_file = '/Users/vahsengrouplaptop/data/phase3/combined_SKB_TPC_ntuples/%s_%s_LUMI_study_new.root'%(month, day)
#Following are for DEC studies...generalize later
HER_input = '/Users/vahsengrouplaptop/data/phase3/combined_SKB_TPC_ntuples/Dec_7_HER_study_new.root'
LER_input = '/Users/vahsengrouplaptop/data/phase3/combined_SKB_TPC_ntuples/Dec_7_LER_study_new.root'

module_ids = ["elepaio","tako","palila"]
#module_ids = ["iiwi", "nene", "humu"]
fit_params_LER = {}

    
bin_width = 240

df_HER_input = rp.read_root(HER_input) #gives dataframe to pass into extract_variables
df_LER_input = rp.read_root(LER_input)
df_LUMI_input = rp.read_root(input_file)

df_HER = extract_variables(df_HER_input, module_ids, "HER", bin_width) #gives dataframe that's computed averages over bin_width slices of input df. Variables can be directly passed into combined heuristic for analysis
df_LER = extract_variables(df_LER_input, module_ids, "LER", bin_width)

###Get single beam fit parameters###
mg_LER, legend_LER, fit_params_LER = make_heuristic_plots(df_LER, module_ids)
mg_HER, legend_HER, fit_params_HER = make_heuristic_plots(df_HER, module_ids)

df = extract_variables_LUMI(df_LUMI_input, module_ids, 30)

###Create dataframe of luminosity data avearged over defined bin_width
for module in module_ids:
    df['LER_bg_'+module]=df['I_LER']*df['P_LER']*2.3**2*fit_params_LER[module].Parameter(0)
    df['LER_T_'+module]=df['I_LER']**2*fit_params_LER[module].Parameter(1)/(df['sy_LER']*df['Nb_LER'])
    df['HER_bg_'+module]=df['I_HER']*df['P_HER']*2.3**2*fit_params_HER[module].Parameter(0)
    df['HER_T_'+module]=df['I_HER']**2*fit_params_HER[module].Parameter(1)/(df['sy_HER']*df['Nb_HER'])
    df['single_beam_fit_'+module] = df['LER_bg_'+module] + df['LER_T_'+module] + df['HER_bg_'+module] + df['HER_T_'+module]
    #df['single_beam_fit_'+module] = df['LER_bg_'+module] + df['LER_T_'+module]

def make_inj_plot(df,module):
    df1 = df.loc[df['dlumi']<25]
    df1 = df1.loc[df1['continuous_inj']==1]
    df1 = df1.loc[df1[module+'_mean']>0]
    df1.index = [i for i in range(0,len(df1))]
    ROOT.gStyle.SetOptFit()
    f_lin = ROOT.TF1("f_lin", "[0]+[1]*x" )
    dlumi = array.array('d', df1['dlumi']/10000)
    lumi = array.array('d', df1['lumi']/10000)
    y = array.array('d',(df1[module+'_mean']-df1['single_beam_fit_'+module]))
    #y = array.array('d',(df1[module+'_mean']))
    y_err = array.array('d',df1[module+'_err'])
    gr = ROOT.TGraphErrors(len(df1), lumi, y, dlumi, y_err)
    fit = gr.Fit("f_lin", "S")
    gr.SetTitle("Lumi Cont. Inj. Fit")
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(0.75)
    gr.SetMarkerColor(4)
    return gr

def make_decay_plot(df,module):
    df1 = df.loc[df['decay']==1]
    df1 = df1.loc[df1[module+'_mean']>0]
    df1.index = [i for i in range(0,len(df1))]
    ROOT.gStyle.SetOptFit()
    f_lin = ROOT.TF1("f_lin", "[0]+[1]*x" )
    diff = df1[module+'_mean'] - df1['single_beam_fit_'+module]
    lumi = array.array('d', df1['lumi']/10000)
    dlumi = array.array('d', df1['dlumi']/10000)
    y = array.array('d',diff)
    dy = array.array('d',df1[module+'_err'])
    gr = ROOT.TGraphErrors(len(df1), lumi, y, dlumi, dy)
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(0.75)
    gr.SetMarkerColor(4)
    gr.SetTitle("%s Lumi Decay Fit"%(module))
    fit = gr.Fit("f_lin", "S")
    df['lumi_background'] = fit.Parameter(0)+fit.Parameter(1)*df['lumi']/10000
    return gr

