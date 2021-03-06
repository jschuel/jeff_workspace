import single_beam_analysis as sa
import pandas as pd
import numpy as np
import root_pandas as rp
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import Patch
import matplotlib
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

rc('text', usetex=True)

tpcs = ["kohola", "nene", "iiwi", "honu", "humu", "tako", "elepaio", "palila"]
phis = ["BWD 18", "BWD 90", "BWD 198", "BWD 270", "FWD 22", "FWD 90", "FWD 202", "FWD 270"]
hist, params_HER = sa.create_histograms(11, "HER")
hist, params_LER = sa.create_histograms(12, "LER")
rates_LER = pd.DataFrame.from_dict(params_LER, orient = 'index').T
rates_HER = pd.DataFrame.from_dict(params_HER, orient = 'index').T
rates_LER = rates_LER[['data_rates_BG', 'data_rates_T', 'sim_rates_BG', 'sim_rates_T']]
rates_HER = rates_HER[['data_rates_BG', 'data_rates_T', 'sim_rates_BG', 'sim_rates_T']]
rates_LER.columns = ['LER_BG', 'LER_T', 'LER_MC_BG', 'LER_MC_T']
rates_HER.columns = ['HER_BG', 'HER_T', 'HER_MC_BG', 'HER_MC_T']
rates_LER[['HER_BG', 'HER_T', 'HER_MC_BG', 'HER_MC_T']] = rates_HER[['HER_BG', 'HER_T', 'HER_MC_BG', 'HER_MC_T']]
rates = rates_LER
rates = rates.reindex(tpcs)
#rates = rates[['LER_MC_BG', 'LER_MC_T', 'HER_MC_BG', 'HER_MC_T', 'LER_BG', 'LER_T', 'HER_BG', 'HER_T']]
rates.index = phis
rates = rates[['LER_BG', 'LER_MC_BG', 'LER_T', 'LER_MC_T', 'HER_BG', 'HER_MC_BG', 'HER_T', 'HER_MC_T']]
rates_total = pd.DataFrame()
rates_total['LER total'] = rates[['LER_BG', 'LER_T']].T.sum()
rates_total['MC LER total'] = rates[['LER_MC_BG', 'LER_MC_T']].T.sum()
rates_total['HER total'] = rates[['HER_BG', 'HER_T']].T.sum()
rates_total['MC HER total'] = rates[['HER_MC_BG', 'HER_MC_T']].T.sum()

#patterns =('','/', '','/','','\\','','\\')
patterns = ('/','/','\\','\\','','','','')
patterns_total = ('', '/', '', '\\')

#colors = ['cyan', matplotlib.colors.colorConverter.to_rgba('#00FFFF', alpha=0.5), 'magenta', matplotlib.colors.colorConverter.to_rgba('#FF00FF', alpha=0.5), 'yellow', matplotlib.colors.colorConverter.to_rgba('#FFFF00', alpha = 0.5), 'indigo', matplotlib.colors.colorConverter.to_rgba('#4b0082', alpha = 0.5)]

colors = [matplotlib.colors.colorConverter.to_rgba('#00FFFF', alpha=0.5), matplotlib.colors.colorConverter.to_rgba('#FF00FF', alpha=0.5), matplotlib.colors.colorConverter.to_rgba('#FFFF00', alpha = 0.5), matplotlib.colors.colorConverter.to_rgba('#4b0082', alpha = 0.5), 'cyan', 'magenta', 'yellow', 'indigo']
colors_total = ['red', matplotlib.colors.colorConverter.to_rgba('#FF0000', alpha=0.5), 'blue', matplotlib.colors.colorConverter.to_rgba('#0000FF', alpha = 0.5)]

#legend_elements = [Patch(facecolor=colors[0], edgecolor='black',label = 'LER BeamGas'), Patch(facecolor=colors[2], edgecolor='black',label = 'LER Touschek'), Patch(facecolor=colors[4], edgecolor='black',label = 'HER BeamGas'), Patch(facecolor=colors[6], edgecolor='black',label = 'HER Touschek'), Patch(facecolor=colors[1], edgecolor='black',label = 'LER BeamGas', hatch = '//'), Patch(facecolor=colors[3], edgecolor='black',label = 'LER Touschek', hatch = '\\\\'), Patch(facecolor=colors[5], edgecolor='black',label = 'HER BeamGas', hatch = '//'), Patch(facecolor=colors[7], edgecolor='black',label = 'HER Touschek', hatch = '\\\\')]

legend_elements = [Patch(facecolor=colors[4], edgecolor='black',label = 'LER BeamGas'), Patch(facecolor=colors[5], edgecolor='black',label = 'LER Touschek'), Patch(facecolor=colors[6], edgecolor='black',label = 'HER BeamGas'), Patch(facecolor=colors[7], edgecolor='black',label = 'HER Touschek'), Patch(facecolor=colors[0], edgecolor='black',label = 'LER BeamGas', hatch = '//'), Patch(facecolor=colors[1], edgecolor='black',label = 'LER Touschek', hatch = '\\\\'), Patch(facecolor=colors[2], edgecolor='black',label = 'HER BeamGas', hatch = '//'), Patch(facecolor=colors[3], edgecolor='black',label = 'HER Touschek', hatch = '\\\\')]  

legend_elements_total = [Patch(facecolor=colors_total[0], edgecolor='black',label = 'LER Total'), Patch(facecolor=colors_total[1], edgecolor='black',label = 'MC LER Total', hatch = '//'), Patch(facecolor=colors_total[2], edgecolor='black',label = 'HER Total'), Patch(facecolor=colors_total[3], edgecolor='black',label = 'MC HER Total', hatch = '\\\\')]

def plot(df = rates, patterns = patterns, colors = colors, legend_elements = legend_elements, stacked = False, logy = True):
    ax = plt.figure(figsize=(10, 6)).add_subplot(111)
    if (stacked == True) & (len(df.columns) > 7):
        df = df[['LER_MC_BG', 'LER_MC_T', 'HER_MC_BG', 'HER_MC_T', 'LER_BG', 'LER_T', 'HER_BG', 'HER_T']]
        if logy == True:
            plt.ylim(4e-1,4e1)
        #df['sum'] = df.sum(axis = 1)
        #for col in df.columns:
        #    df[col] = df[col]/df['sum']
        #df = df[['LER_MC_BG', 'LER_MC_T', 'HER_MC_BG', 'HER_MC_T', 'LER_BG', 'LER_T', 'HER_BG', 'HER_T']]
    df.plot(kind = 'bar', stacked = stacked, color = colors, linewidth = 1, logy = logy, edgecolor = 'black', ax = ax, legend = False, width = 0.75)
    bars = ax.patches
    hatches = [p for p in patterns for i in range(len(df))]
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
    ax.legend(handles = legend_elements, title = r'Data $\qquad \qquad \qquad \qquad \qquad \quad$MC', ncol = 2)
    plt.tight_layout()
    plt.rc('legend', fontsize=12)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    plt.rc('axes', labelsize=18)
    plt.rc('axes', titlesize=18)
    plt.ylabel("Rate [Hz]")
    plt.axes().yaxis.set_minor_locator(AutoMinorLocator())
    matplotlib.rc('xtick', labelsize=20)
    matplotlib.rc('ytick', labelsize=20)

plot(df = rates_total, patterns = patterns_total, colors = colors_total, legend_elements = legend_elements_total, stacked = True, logy = False)
