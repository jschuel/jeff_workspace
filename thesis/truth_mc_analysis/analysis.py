import ROOT
import root_pandas as rp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import array
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.collections import PatchCollection

class lepton_analysis:
    def __init__(self):
        self.data = self.get_data_with_cuts()
        self.signal = self.get_signal()
        self.bg = self.get_background()
    
    def load_data(self, input_file = '~/workspace/thesis_data/mc_with_deltaT.root', key = 'truth_lepton'):
        data = rp.read_root(input_file, key=key)
        return data

    def get_data_with_cuts(self, input_file = '~/workspace/thesis_data/mc_with_deltaT.root', key = 'truth_lepton'):
        data = self.load_data()
        data_post_cut = data.loc[data['Mbc'] < 5.245]
        return data_post_cut

    def get_signal(self):
        cut_data = self.get_data_with_cuts()
        signal = cut_data.loc[np.abs(cut_data['genMotherPDG']) == 511]
        return signal

    def get_background(self):
        cut_data = self.get_data_with_cuts()
        background = cut_data.loc[np.abs(cut_data['genMotherPDG']) != 511]
        return background

    def process_deltaT_and_flavor_tag(self, data): #only use with daughter class
        data['ncands'] = [len(data.loc[data['__event__'] == event]) for event in data['__event__']]
        data = data.loc[data['ncands']==2]
        data['PDG_sum'] = [np.abs(data.loc[data['__event__']==event]['mcPDG'].sum()) for event in data['__event__']] #sums PDG code for flavor tagging
        data['OF'] = 0 #creates opposite flavor_tag column
        data['SF'] = 0 #creates same flavor_tag columns
        OF_index = data.loc[data['PDG_sum']<10].index.to_numpy() #gets indices where events are opposite flavor
        SF_index = data.loc[data['PDG_sum']>10].index.to_numpy() # gets indices where events are same flavor
        data['OF'][OF_index] = 1 #tags all opposite flavor events with 1
        data['SF'][SF_index] = 1
        data['sumZ'] = [np.abs(data.loc[data['__event__']==event]['z'].sum()) for event in data['__event__']]
        data['deltaZ'] = [np.abs(data.loc[data['__event__'] == event]['z'].diff().iloc[1]) for event in data['__event__']]
        data['sumT'] = data['sumZ']/(0.425*3e10)*1e12
        data['deltaT'] = data['deltaZ']/(0.425*3e10)*1e12
        data = data.drop(columns = ['ncands', 'PDG_sum']) #drops filler columns used to flavor tag
        return data

    def bin_data_and_compute_asymmetry(self, data):
        data = data.loc[data['deltaT'].duplicated() == False]
        data_reduced = data[['OF', 'SF', 'deltaT', 'sumT']]
        data_reduced['bin_num'] = 0
        #bin according to Go's convention for now
        data_reduced['bin_num'][data_reduced.loc[(data_reduced['deltaT']>=0.5) & (data_reduced['deltaT']< 1)].index.to_numpy()]=1 #Bin_Num according to Go's convention
        for i in range(2, 8):
            data_reduced['bin_num'][data_reduced.loc[(data_reduced['deltaT']>=i-1) & (data_reduced['deltaT']< i)].index.to_numpy()]=i
        data_reduced['bin_num'][data_reduced.loc[(data_reduced['deltaT']>=7) & (data_reduced['deltaT']< 9)].index.to_numpy()]=8
        data_reduced['bin_num'][data_reduced.loc[(data_reduced['deltaT']>=9) & (data_reduced['deltaT']< 13)].index.to_numpy()]=9
        data_reduced['bin_num'][data_reduced.loc[(data_reduced['deltaT']>=13) & (data_reduced['deltaT']< 20)].index.to_numpy()]=10
        data_reduced['asymmetry'] = 0
        asymmetry = []
        for i in range(0,11):
            asymmetry.append((data_reduced.loc[data_reduced['bin_num']==i]['OF'].sum()-data_reduced.loc[data_reduced['bin_num']==i]['SF'].sum())/(data_reduced.loc[data_reduced['bin_num']==i]['OF'].sum()+data_reduced.loc[data_reduced['bin_num']==i]['SF'].sum()))
            index = data_reduced.loc[data_reduced['bin_num'] == i].index.to_numpy()
            data_reduced['asymmetry'][index] = asymmetry[i]
        return data_reduced

    def get_data_for_plotting(self, data):
        initial_data = self.bin_data_and_compute_asymmetry(data)
        d = initial_data.loc[initial_data['bin_num'].duplicated() == False] #get unique asymmetries
        d = d.sort_values(by = 'bin_num')
        d['x'] = np.array([0.25, 0.75, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 8, 11, 16.5])
        d['xerr'] = np.array([0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 2, 3.5])
        d['y'] = d['asymmetry'].to_numpy()
        d['yerr'] = yerr = np.array([0 for i in range(0,len(d))])
        d['y_go'] = np.array([1.013, 0.916, 0.699, 0.339, -0.136, -0.634, -0.961, -0.974, -0.675, 0.089, 0.243]) #Go's values
        d['yerr_go'] = np.array([0.028, 0.022, 0.038, 0.056, 0.075, 0.084, 0.077, 0.080, 0.109, 0.193, 0.435]) #Go's values
        d = d.drop(columns = ['bin_num', 'deltaT', 'asymmetry'])
        return d

class truth_B_analysis:
    def __init__(self):
        self.data = self.load_data()
    
    def load_data(self, input_file = '~/workspace/thesis_data/mc_with_deltaT.root', key = 'truthB_lab'):
        data = rp.read_root(input_file, key=key)
        return data

    def process_deltaT_and_flavor_tag(self, data): #only use with daughter class
        data['ncands'] = [len(data.loc[data['__event__'] == event]) for event in data['__event__']]
        data = data.loc[data['ncands']==2]
        data['PDG_sum'] = [np.abs(data.loc[data['__event__']==event]['mcPDG'].sum()) for event in data['__event__']] #sums PDG code for flavor tagging
        data['OF'] = 0 #creates opposite flavor_tag column
        data['SF'] = 0 #creates same flavor_tag columns
        OF_index = data.loc[data['PDG_sum']==0].index.to_numpy() #gets indices where events are opposite flavor
        SF_index = data.loc[data['PDG_sum']!=0].index.to_numpy() # gets indices where events are same flavor
        data['OF'][OF_index] = 1 #tags all opposite flavor events with 1
        data['SF'][SF_index] = 1
        data['sumZ'] = [np.abs(data.loc[data['__event__']==event]['mcZ'].sum()) for event in data['__event__']]
        data['deltaZ'] = [np.abs(data.loc[data['__event__'] == event]['mcZ'].diff().iloc[1]) for event in data['__event__']]
        data['sumT'] = [np.abs(data.loc[data['__event__']==event]['mcZ'].sum()) for event in data['__event__']]
        data['deltaT'] = data['deltaZ']/(0.425*3e10)*1e12
        data = data.drop(columns = ['ncands', 'PDG_sum']) #drops filler columns used to flavor tag
        return data

    def bin_data_and_compute_asymmetry(self, data, col = 'deltaT'):
        data = data.loc[data[col].duplicated() == False]
        data_reduced = data[['OF', 'SF', col, 'sumT']]
        data_reduced['bin_num'] = 0
        #bin according to Go's convention for now
        data_reduced['bin_num'][data_reduced.loc[(data_reduced[col]>=0.5) & (data_reduced[col]< 1)].index.to_numpy()]=1 #Bin_Num according to Go's convention
        for i in range(2, 8):
            data_reduced['bin_num'][data_reduced.loc[(data_reduced[col]>=i-1) & (data_reduced[col]< i)].index.to_numpy()]=i
        data_reduced['bin_num'][data_reduced.loc[(data_reduced[col]>=7) & (data_reduced[col]< 9)].index.to_numpy()]=8
        data_reduced['bin_num'][data_reduced.loc[(data_reduced[col]>=9) & (data_reduced[col]< 13)].index.to_numpy()]=9
        data_reduced['bin_num'][data_reduced.loc[(data_reduced[col]>=13) & (data_reduced[col]< 20)].index.to_numpy()]=10
        data_reduced['asymmetry'] = 0
        asymmetry = []
        for i in range(0,11):
            asymmetry.append((data_reduced.loc[data_reduced['bin_num']==i]['OF'].sum()-data_reduced.loc[data_reduced['bin_num']==i]['SF'].sum())/(data_reduced.loc[data_reduced['bin_num']==i]['OF'].sum()+data_reduced.loc[data_reduced['bin_num']==i]['SF'].sum()))
            index = data_reduced.loc[data_reduced['bin_num'] == i].index.to_numpy()
            data_reduced['asymmetry'][index] = asymmetry[i]
        return data_reduced

    def get_data_for_plotting(self, data, col = 'deltaT'):
        initial_data = self.bin_data_and_compute_asymmetry(data, col)
        d = initial_data.loc[initial_data['bin_num'].duplicated() == False] #get unique asymmetries
        d = d.sort_values(by = 'bin_num')
        d['x'] = np.array([0.25, 0.75, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 8, 11, 16.5])
        d['xerr'] = np.array([0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 2, 3.5])
        d['y'] = d['asymmetry'].to_numpy()
        d['yerr'] = yerr = np.array([0 for i in range(0,len(d))])
        d['y_go'] = np.array([1.013, 0.916, 0.699, 0.339, -0.136, -0.634, -0.961, -0.974, -0.675, 0.089, 0.243]) #Go's values
        d['yerr_go'] = np.array([0.028, 0.022, 0.038, 0.056, 0.075, 0.084, 0.077, 0.080, 0.109, 0.193, 0.435]) #Go's values
        d = d.drop(columns = ['bin_num', col, 'asymmetry'])
        return d

class process_lepton_data(lepton_analysis): #only call to add deltaT branches to root file. Takes a while to run
    def __init__(self, input_data = '~/workspace/thesis_data/mc_with_reco.root', output_data = '~/workspace/thesis_data/mc_with_deltaT.root', key = 'truth_lepton'):
        super().__init__()
        self.data = self.process_deltaT_and_flavor_tag(self.get_data_with_cuts(input_file=input_data, key=key))
        self.data.to_root(output_data, key = key)

class process_B_data(truth_B_analysis): #only call to add deltaT branches to root file. Takes a while to run
    def __init__(self, input_data = '~/workspace/thesis_data/mc_with_reco.root', output_data = '~/workspace/thesis_data/mc_with_deltaT.root', key = 'truthB_lab'):
        super().__init__()
        self.data = self.process_deltaT_and_flavor_tag(self.load_data(input_file=input_data, key=key))
        self.data.to_root(output_data, key = key, mode = 'a')

def plot_data():
    la = lepton_analysis() #initialize lepton analysis class
    lepton_data = la.get_data_with_cuts()
    lepton_signal = la.get_signal()
    B = truth_B_analysis() #initialize truth B analysis class
    B_data = B.load_data()
    lepton_plot_df = la.get_data_for_plotting(lepton_data) #la.data is the data pulled from initializing the la class
    lepton_plot_signal = la.get_data_for_plotting(lepton_signal)
    B_plot_df_reco = B.get_data_for_plotting(B_data, col = 'deltaT') #uses deltaT reconstructed from truth B decay vertices
    B_plot_df = B.get_data_for_plotting(B_data, col = 'MCDeltaT') #uses built in basf2 deltaT function
    plt.errorbar(lepton_plot_df['x'], lepton_plot_df['y'], lepton_plot_df['yerr'], lepton_plot_df['xerr'], 'o', markersize = 0, label = 'delta from reconstructed truth lepton production vertices with Mbc<5.245')
    plt.errorbar(lepton_plot_signal['x'], lepton_plot_signal['y'], lepton_plot_signal['yerr'], lepton_plot_signal['xerr'], 'o', markersize = 0, label = 'Signal only deltaT computed from reconstructed truth lepton production vertices with Mbc<5.245')
    plt.errorbar(B_plot_df_reco['x'], B_plot_df_reco['y'], B_plot_df_reco['yerr'], B_plot_df_reco['xerr'], 'o', markersize = 0, label = 'deltaT computed from truth B decay vertices')
    plt.errorbar(B_plot_df['x'], B_plot_df['y'], B_plot_df['yerr'], B_plot_df['xerr'], 'o', markersize = 0, label = 'deltaT computed from built in basf2 truth deltaT')
    plt.errorbar(lepton_plot_df['x'], lepton_plot_df['y_go'], lepton_plot_df['yerr_go'], lepton_plot_df['xerr'], 'o', markersize = 0, label = "Go's data")
    plt.xlabel('deltaT [ps]')
    plt.ylabel('Asymmetry(deltaT)')
    plt.legend()
    plt.show()

def plot_with_fits():
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    la = lepton_analysis() #initialize lepton analysis class
    lepton_signal = la.get_signal()
    B = truth_B_analysis() #initialize truth B analysis class
    B_data = B.load_data()
    lepton_plot_df = la.get_data_for_plotting(lepton_signal) #lepton signal daata
    B_plot_df = B.get_data_for_plotting(B_data, col = 'deltaT') #uses deltaT reconstructed from truth B decay vertices
    lepton_fits, B_fits = fit_all_cases()
    lepton_qm_max = lepton_fits['QM'] + 20*lepton_fits['QM_err']
    B_qm_max = B_fits['QM'] + 20*B_fits['QM_err']
    lepton_qm_min = lepton_fits['QM'] - 20*lepton_fits['QM_err']
    B_qm_min = B_fits['QM'] - 20*B_fits['QM_err']
    lepton_sd_max = lepton_fits['SD'] + 20*lepton_fits['SD_err']
    B_sd_max = B_fits['SD'] + 20*B_fits['SD_err']
    lepton_sd_min = lepton_fits['SD'] - 20*lepton_fits['SD_err']
    B_sd_min = B_fits['SD'] - 20*B_fits['SD_err']
    #x_plt = np.linspace(0,20, 201)
    y_lepton_qm_min = np.cos(lepton_qm_min*lepton_plot_df['x'])
    y_lepton_qm_max = np.cos(lepton_qm_max*lepton_plot_df['x'])
    y_lepton_sd_min = np.cos(lepton_sd_min*lepton_plot_df['x'])
    y_lepton_sd_max = np.cos(lepton_sd_max*lepton_plot_df['x'])
    y_B_qm_min = np.cos(B_qm_min*B_plot_df['x'])
    y_B_qm_max = np.cos(B_qm_max*B_plot_df['x'])
    y_B_sd_min = np.cos(B_sd_min*B_plot_df['x'])
    y_B_sd_max = np.cos(B_sd_max*B_plot_df['x'])
    #y_fit = np.cos(fit['truth_MC']*x_data)
    #residuals = y_data - y_fit
    fig, ax = plt.subplots(nrows = 2, ncols = 2)
    ax[0,0].errorbar(lepton_plot_df['x'], lepton_plot_df['y'], lepton_plot_df['yerr'], lepton_plot_df['xerr'], 'o', label = r'Data from reconstructed truth lepton prod. vertices with $M_{bc}<5.245$')
    make_error_boxes(ax[0,0], lepton_plot_df['x'], y_lepton_qm_min, lepton_plot_df['xerr'], y_lepton_qm_max)
    ax[0,0].add_patch(Rectangle((lepton_plot_df['x'][0] - lepton_plot_df['xerr'][0], y_lepton_qm_min[0]), 2*lepton_plot_df['xerr'][0], y_lepton_qm_max[0]-y_lepton_qm_min[0], edgecolor = 'green', facecolor = 'yellow', alpha = 0.5, label = r'1$\sigma$ error region for QM model fit. $\Delta m_d =$ %s $\pm$ %s ps$^{-1}$'%(float('%.3g' % lepton_fits['QM']), float('%.3g' % lepton_fits['QM_err']))))
    ax[0,0].legend()
    ax[0,0].set_title('QM model fit to truth lepton prod. vertex results')
    ax[0,0].set_ylabel(r"$\mathcal{A}(\Delta t)$")
    ax[0,1].errorbar(lepton_plot_df['x'], lepton_plot_df['y'], lepton_plot_df['yerr'], lepton_plot_df['xerr'], 'o', label = r'Data from reconstructed truth lepton prod. vertices with $M_{bc}<5.245$')
    make_error_boxes(ax[0,1], lepton_plot_df['x'], y_lepton_sd_min, lepton_plot_df['xerr'], y_lepton_sd_max)
    ax[0,1].add_patch(Rectangle((lepton_plot_df['x'][0] - lepton_plot_df['xerr'][0], y_lepton_sd_min[0]), 2*lepton_plot_df['xerr'][0], y_lepton_sd_max[0]-y_lepton_sd_min[0], edgecolor = 'green', facecolor = 'yellow', alpha = 0.5, label = r'1$\sigma$ error region for SD model fit. $\Delta m_d =$ %s $\pm$ %s ps$^{-1}$'%(float('%.3g' % lepton_fits['SD']), float('%.3g' % lepton_fits['SD_err']))))
    ax[0,1].legend()
    ax[0,1].set_title('SD model fit to truth lepton prod. vertex results')
    ax[0,1].set_ylabel(r"$\mathcal{A}(\Delta t)$")
    ax[1,0].errorbar(B_plot_df['x'], B_plot_df['y'], B_plot_df['yerr'], B_plot_df['xerr'], 'o', label = 'Data from reconstructed truth B decay vertices')
    make_error_boxes(ax[1,0], B_plot_df['x'], y_B_qm_min, B_plot_df['xerr'], y_B_qm_max)
    ax[1,0].add_patch(Rectangle((B_plot_df['x'][0] - B_plot_df['xerr'][0], y_B_qm_min[0]), 2*B_plot_df['xerr'][0], y_B_qm_max[0]-y_B_qm_min[0], edgecolor = 'green', facecolor = 'yellow', alpha = 0.5, label = r'1$\sigma$ error region for QM model fit. $\Delta m_d =$ %s $\pm$ %s ps$^{-1}$'%(float('%.3g' % B_fits['QM']), float('%.3g' % B_fits['QM_err']))))
    ax[1,0].legend()
    ax[1,0].set_title('QM model fit to truth B decay vertex results')
    ax[1,0].set_ylabel(r"$\mathcal{A}(\Delta t)$")
    ax[1,1].errorbar(B_plot_df['x'], B_plot_df['y'], B_plot_df['yerr'], B_plot_df['xerr'], 'o', label = 'Data from reconstructed truth B decay vertices')
    make_error_boxes(ax[1,1], B_plot_df['x'], y_B_sd_min, B_plot_df['xerr'], y_B_sd_max)
    ax[1,1].add_patch(Rectangle((B_plot_df['x'][0] - B_plot_df['xerr'][0], y_B_sd_min[0]), 2*B_plot_df['xerr'][0], y_B_sd_max[0]-y_B_sd_min[0], edgecolor = 'green', facecolor = 'yellow', alpha = 0.5, label = r'1$\sigma$ error region for SD model fit. $\Delta m_d =$ %s $\pm$ %s ps$^{-1}$'%(float('%.3g' % B_fits['SD']), float('%.3g' % B_fits['SD_err']))))
    ax[1,1].legend()
    ax[1,1].set_title('SD model fit to truth B decay vertex results')
    ax[1,1].set_ylabel(r"$\mathcal{A}(\Delta t)$")
    fig.show()
    #ax[1].bar(x_data,residuals,width=2*x_err_data, fill=False)
    #ax[1].axhline(y=0, xmin = 0, xmax = 20)
    #ax[1].set_ylabel(r"$\mathcal{A}_{Truth_{MC}}-\mathcal{A}_{Fit}$")
    #ax[1].set_xlabel(r"$\Delta t$")

def root_fit(plot_data, col = 'deltaT'): #fits cosine curve to data (asymmetry computed from binning)
    plot_data = plot_data.loc[plot_data[col]<=20]
    plot_data = plot_data.loc[plot_data['sumT']<=250]
    x1 = array.array('d', plot_data[col])
    x2 = array.array('d', plot_data['sumT'])
    y = array.array('d', plot_data['asymmetry'])
    gr = ROOT.TGraph(len(x1), x1, y)
    f = ROOT.TF1('f1', 'cos([0]*x)', 0, 20)
    f.SetParLimits(0,0.1,0.9)
    gr.Fit('f1', 'S')
    f2 = ROOT.TF2("f2","0.5*(cos([0]*x)+cos([0]*y))", 0, 20, 0, 100)
    f2.SetParLimits(0,0.1,0.9)
    gr2 = ROOT.TGraph2D(len(x1), x1, x2, y)
    gr2.Fit('f2', 'S')
    fitdict = {"QM":f.GetParameter(0), "QM_err":f.GetParError(0), "SD":f2.GetParameter(0), "SD_err":f.GetParError(0)}
    return fitdict
    
def fit_all_cases(): #gives fit results for each MC test sample
    la = lepton_analysis()
    B = truth_B_analysis()
    l_data = la.get_data_with_cuts()
    B_data = B.load_data()
    l_plot_data = la.bin_data_and_compute_asymmetry(l_data)
    B_plot_data = B.bin_data_and_compute_asymmetry(B_data)
    lepton_fits = root_fit(l_plot_data)
    B_fits = root_fit(B_plot_data)
    return lepton_fits, B_fits

def make_error_boxes(ax, xdata, ymin, xerror, ymax, facecolor='yellow',
                     edgecolor='green', alpha=0.5):
    # Loop over data points; create box from errors at each point
    errorboxes = [Rectangle((x - xe, y), 2*xe, ye-y) for x, y, xe, ye in zip(xdata, ymin, xerror, ymax)]
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, facecolor=facecolor, alpha=alpha,
                         edgecolor=edgecolor)
    # Add collection to axes
    ax.add_collection(pc)

plot_with_fits()
