import pandas as pd
import root_pandas as rp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import ROOT
from ROOT import TVector3
from array import array
import matplotlib
from matplotlib.lines import Line2D
from matplotlib import rc
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
rc('text', usetex=False)
pd.set_option('mode.chained_assignment', None) #remove copy warning
plt.rc('legend', fontsize=12)
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)
plt.rc('axes', labelsize=16)
plt.rc('axes', titlesize=16)

class analysis:
    def __init__(self, E_cut = {'palila': 8, 'iiwi': 8.8, 'tako': 6, 'nene': 5.5, 'elepaio': 6, 'humu': 6.6}, skb_file= "/home/jeff/data/phase3/spring_2020/05-09-20/combined_ntuples/05-09_whole_study_2021_update.root", tpc_path =  '/home/jeff/data/phase3/spring_2020/05-09-20/tpc_root_files/', recoils_only = True, fei4_restrict = True, data_machine_params = {'I_HER' : 510, 'I_LER' : 510, 'P_LER': 3e-8, 'P_HER': 1.4e-8, 'sy_LER' : 60, 'sy_HER' : 35, 'nb_LER' : 783, 'nb_HER' : 783, 'lumi' : 1.1}, MC_machine_params = {'I_HER' : 1000, 'I_LER' : 1200, 'sy_LER' : 37, 'sy_HER' : 36, 'nb_LER' : 1576, 'nb_HER' : 1576, 'lumi' : 25}, bin_width = 120, fit_package = 'scipy'): # for fit use 'scipy' or 'root'   Old skb_file: 05-09_whole_study_even_newerest.root
        self.E_cut = E_cut
        self.tpcs = {'z = -14m': 'elepaio', 'z = -8.0m': 'tako', 'z = -5.6m': 'palila', 'z = +6.6m': 'iiwi', 'z = +14m': 'nene', 'z = +16m': 'humu'}
        self.skb_file = skb_file
        self.tpc_path = tpc_path
        self.recoils_only = recoils_only
        self.fei4_restrict = fei4_restrict
        self.fit_package = fit_package
        self.raw_tpc_data = self.get_tpc_data()
        self.raw_study_data = self.get_raw_study_data()
        self.data = self.combine_tpc_and_SKB_data()
        self.binned_data = self.bin_data(nsecs = bin_width)
        self.SB_fit_params = self.generate_fits()
        self.lumi_fit_params = self.fit_lumi()[0]
        self.MC_bgTypes = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all', 'RBB_Lumi', 'twoPhoton_Lumi']
        self.MC_machine_params = MC_machine_params
        self.MC_data = self.get_MC_data()
        self.MC_rates = self.get_MC_rates()
        self.data_machine_params = data_machine_params

    def get_tpc_data(self):
        data = {}
        for tpc in self.tpcs.values():
            data[tpc] = rp.read_root(self.tpc_path + "%s_all_newester8.root"%(tpc))
            if isinstance(self.E_cut,dict):
                data[tpc] = data[tpc].loc[data[tpc]['track_energy']>=(self.E_cut[tpc])]
            else:
                data[tpc] = data[tpc].loc[data[tpc]['track_energy']>=self.E_cut]
            if self.recoils_only:
                data[tpc] = data[tpc].loc[data[tpc]['is_recoil'] == 1]
            data[tpc]['ts'] = data[tpc]['timestamp_start'].astype('int')
        return data

    def get_raw_study_data(self):
        study_data = rp.read_root(self.skb_file)
        study_data['ts'] = study_data['ts'].astype('int')
        tpc_data = self.raw_tpc_data
        dfs = {}
        for tpc in self.tpcs.values():
            if isinstance(self.E_cut,dict):
                dfs[tpc] = tpc_data[tpc].loc[tpc_data[tpc]['track_energy']>=(self.E_cut[tpc])]
            else:
                dfs[tpc] = tpc_data[tpc].loc[tpc_data[tpc]['track_energy']>=(self.E_cut)]
            study_data[tpc+'_neutrons'] = 0
            pv = dfs[tpc].pivot_table(index=['ts'],aggfunc='size')
            pv = pv.loc[pv.index.isin(study_data['ts'])]
            counts = pd.DataFrame()
            counts['ts'] = pv.index
            counts['rate'] = pv.array
            index = study_data.loc[study_data['ts'].astype('int').isin(counts['ts'])].index.to_numpy()
            study_data[tpc+'_neutrons'][index] = counts['rate'].to_numpy()
        return study_data

    def combine_tpc_and_SKB_data(self):
        skb = self.raw_study_data
        tpc_data = self.raw_tpc_data
        tpc_reduced = {} # tpc data within study fills
        for tpc in self.tpcs.values():
            tpc_data[tpc]['count'] = 1
            tpc_reduced[tpc] = tpc_data[tpc].loc[(tpc_data[tpc]['ts'].isin(skb['ts']))]
            tpc_reduced[tpc].index = [i for i in range(0,len(tpc_reduced[tpc]))]
            skb['%s_neutrons'%(tpc)] = 0 #reset these branches to fill with entrys from tpc_reduced
            index = skb.loc[skb['ts'].isin(tpc_reduced[tpc]['ts'])].index.to_numpy()
            skb['%s_neutrons'%(tpc)][index] = tpc_reduced[tpc].groupby(tpc_reduced[tpc]['ts']).sum()['count'].to_numpy()
            tmp = skb.loc[skb['ts'].isin(tpc_reduced[tpc]['ts'])]
            I_LER = []
            I_HER = []
            P_LER = []
            P_HER = []
            Sy_LER = []
            Sy_HER = []
            Sx_LER = []
            Sx_HER = []
            Nb_LER = []
            Nb_HER = []
            Lumi = []
            LER_flag = []
            HER_flag = []
            Lumi_flag = []
            cinj_flag = []
            decay_flag = []
            for Il, Ih, Pl, Ph, Syl, Syh, Sxl, Sxh, Nbl, Nbh, lumi, rate, fler, fher, flum, fci, fd in zip(tmp['I_LER'], tmp['I_HER'], tmp['P_LER'], tmp['P_HER'], tmp['Sy_LER'], tmp['Sy_HER'], tmp['Sx_LER'], tmp['Sx_HER'], tmp['Nb_LER'], tmp['Nb_HER'], tmp['ECL_lumi'], tmp['%s_neutrons'%(tpc)], tmp['LER_study_flag'], tmp['HER_study_flag'], tmp['Lumi_study_flag'],tmp['Cont_inj_flag'], tmp['Decay_flag']):
                if rate > 1:
                    for i in range(0,rate):
                        I_LER.append(Il)
                        I_HER.append(Ih)
                        P_LER.append(Pl)
                        P_HER.append(Ph)
                        Sy_LER.append(Syl)
                        Sy_HER.append(Syh)
                        Sx_LER.append(Sxl)
                        Sx_HER.append(Sxh)
                        Nb_LER.append(Nbl)
                        Nb_HER.append(Nbh)
                        Lumi.append(lumi)
                        LER_flag.append(fler)
                        HER_flag.append(fher)
                        Lumi_flag.append(flum)
                        cinj_flag.append(fci)
                        decay_flag.append(fd)
                else:
                    I_LER.append(Il)
                    I_HER.append(Ih)
                    P_LER.append(Pl)
                    P_HER.append(Ph)
                    Sy_LER.append(Syl)
                    Sy_HER.append(Syh)
                    Sx_LER.append(Sxl)
                    Sx_HER.append(Sxh)
                    Nb_LER.append(Nbl)
                    Nb_HER.append(Nbh)
                    Lumi.append(lumi)
                    LER_flag.append(fler)
                    HER_flag.append(fher)
                    Lumi_flag.append(flum)
                    cinj_flag.append(fci)
                    decay_flag.append(fd)
            tpc_reduced[tpc]['I_LER'] = I_LER
            tpc_reduced[tpc]['I_HER'] = I_HER
            tpc_reduced[tpc]['P_LER'] = P_LER
            tpc_reduced[tpc]['P_HER'] = P_HER
            tpc_reduced[tpc]['Sy_LER'] = Sy_LER
            tpc_reduced[tpc]['Sy_HER'] = Sy_HER
            tpc_reduced[tpc]['Sx_LER'] = Sx_LER
            tpc_reduced[tpc]['Sx_HER'] = Sx_HER
            tpc_reduced[tpc]['Nb_LER'] = Nb_LER
            tpc_reduced[tpc]['Nb_HER'] = Nb_HER
            tpc_reduced[tpc]['Lumi'] = Lumi
            tpc_reduced[tpc]['LER_flag'] = LER_flag
            tpc_reduced[tpc]['HER_flag'] = HER_flag
            tpc_reduced[tpc]['Lumi_flag'] = Lumi_flag
            tpc_reduced[tpc]['Cont_inj_flag'] = cinj_flag
            tpc_reduced[tpc]['Decay_flag'] = decay_flag

            LER_index = tpc_reduced[tpc].loc[(tpc_reduced[tpc]['LER_flag'] == 1) &
                                             ((tpc_reduced[tpc]['Cont_inj_flag'] == 1) |
                                              (tpc_reduced[tpc]['Decay_flag'] == 1))].index.to_numpy()
            tpc_reduced[tpc]['LER_flag'] = 0
            tpc_reduced[tpc]['LER_flag'][LER_index]=1 #reindex study flag to only be good continuous injection and decay periods
            
            HER_index = tpc_reduced[tpc].loc[(tpc_reduced[tpc]['HER_flag'] == 1) &
                                             ((tpc_reduced[tpc]['Cont_inj_flag'] == 1) |
                                              (tpc_reduced[tpc]['Decay_flag'] == 1))].index.to_numpy()
            tpc_reduced[tpc]['HER_flag'] = 0
            tpc_reduced[tpc]['HER_flag'][HER_index]=1 #reindex study flag to only be good continuous injection and decay periods

            Lumi_index = tpc_reduced[tpc].loc[(tpc_reduced[tpc]['Lumi_flag'] == 1) &
                                             ((tpc_reduced[tpc]['Cont_inj_flag'] == 1) |
                                              (tpc_reduced[tpc]['Decay_flag'] == 1))].index.to_numpy()
            tpc_reduced[tpc]['Lumi_flag'] = 0
            tpc_reduced[tpc]['Lumi_flag'][Lumi_index]=1 #reindex study flag to only be good continuous injection and decay periods
            
        return tpc_reduced

    def bin_data(self,nsecs): #nsecs is number of seconds per bin
        data = self.data
        avgs = {}
        errs = {}
        for tpc in self.tpcs.values():
            avgs[tpc] = data[tpc].groupby(pd.cut(data[tpc]['ts'], bins = np.linspace(data[tpc]['ts'].min(),data[tpc]['ts'].max(),int((data[tpc]['ts'].max()-data[tpc]['ts'].min())/nsecs)))).mean() #bin into bins of width nsecs
            errs[tpc] = data[tpc].groupby(pd.cut(data[tpc]['ts'], bins = np.linspace(data[tpc]['ts'].min(),data[tpc]['ts'].max(),int((data[tpc]['ts'].max()-data[tpc]['ts'].min())/nsecs)))).std() #bin into bins of width nsecs
            for col in ['ts', 'I_LER', 'I_HER', 'P_LER', 'P_HER', 'Sy_LER', 'Sy_HER', 'Sx_LER', 'Sx_HER', 'Nb_LER', 'Nb_HER', 'Lumi']:
                avgs[tpc][col+'_err'] = errs[tpc][col]
            avgs[tpc]['rate'] = data[tpc].groupby(pd.cut(data[tpc]['ts'], bins = np.linspace(data[tpc]['ts'].min(),data[tpc]['ts'].max(),int((data[tpc]['ts'].max()-data[tpc]['ts'].min())/nsecs))))['count'].sum()/nsecs
            avgs[tpc]['rate_err'] = np.sqrt(data[tpc].groupby(pd.cut(data[tpc]['ts'], bins = np.linspace(data[tpc]['ts'].min(),data[tpc]['ts'].max(),int((data[tpc]['ts'].max()-data[tpc]['ts'].min())/nsecs))))['count'].sum())/nsecs
            if self.fit_package == 'scipy':
                avgs[tpc] = avgs[tpc].dropna() #drop rows without entries
        return avgs

    def generate_fits(self):
        data = self.binned_data
        if self.fit_package == 'scipy':
            def fit(tpc,y,yerr,x1,x1err,x2,x2err):
                def func(data,a,b,c):
                    y = a*data[:,0]-b+c*data[:,1]
                    return y
                data = pd.DataFrame()
                data['x1'] = x1
                data['x2'] = x2
                data['y'] = y
                data['yerr'] = yerr
                fit_data = data[['x1','x2']].to_numpy()
                ydata = data['y'].to_numpy()
                yerror = data['yerr'].to_numpy()
                print(ydata)
                print(yerror)
                if len(ydata)>3:
                    try:
                        params, pcov = optimize.curve_fit(func, fit_data, ydata, sigma=yerror, absolute_sigma = True, bounds = ([0,0,0], [1e5,1e-2,10]),maxfev=100000)
                        B0 = params[0]
                        B1 = params[1]
                        T = params[2]
                        B0_err = np.sqrt(np.diag(pcov))[0]
                        B1_err = np.sqrt(np.diag(pcov))[1]
                        T_err = np.sqrt(np.diag(pcov))[2]
                    except ValueError:
                        B0 = 0
                        B1 = 0
                        T = 0
                        B0_err = 1
                        B1_err = 1
                        T_err = 1
                else:
                    B0 = 0
                    B1 = 0
                    T = 0
                    B0_err = 1
                    B1_err = 1
                    T_err = 1
                return B0, B1, T, B0_err, B1_err, T_err
        else:
            def fit(tpc,y,yerr,x1,x1err,x2,x2err):
                f = ROOT.TF2("f2","[0]*x - [1] + [2]*y")
                f.SetParLimits(0, 0, 1e5)
                f.SetParLimits(1,0,1e-2)
                f.SetParLimits(2,0,10)
                gr = ROOT.TGraph2DErrors(len(x1), x1, x2, y, x1err, x2err, yerr)
                gr.Fit(f, 'SEM')
                B0 = f.GetParameter(0)
                B1 = f.GetParameter(1)
                T = f.GetParameter(2)
                B0_err = f.GetParError(0)
                B1_err = f.GetParError(1)
                T_err = f.GetParError(2)
                return B0, B1, T, B0_err, B1_err, T_err
        
        '''
        def fit(tpc,y,yerr,x1,x1err,x2,x2err):
            f = ROOT.TF2("f2","[0]*x + [1]*y")
            f.SetParLimits(0,0,10)
            gr = ROOT.TGraph2DErrors(len(x1), x1, x2, y, x1err, x2err, yerr)
            gr.Fit(f, 'SEM')
            B1 = f.GetParameter(0)
            T = f.GetParameter(1)
            B1_err = f.GetParError(0)
            T_err = f.GetParError(1)
            return B1, T, B1_err, T_err 
        
        def fit(tpc,y,yerr,x1,x1err):
            f = ROOT.TF1("f","[0] + [1]*x")
            f.SetParLimits(0,0,10)
            f.SetParLimits(1,0,10)
            gr = ROOT.TGraphErrors(len(x1), x1,  y, x1err, yerr)
            gr.Fit(f, 'SEM')
            B1 = f.GetParameter(0)
            T = f.GetParameter(1)
            B1_err = f.GetParError(0)
            T_err = f.GetParError(1)
            return B1, T, B1_err, T_err
        '''
        
        fits = {}
        for tpc in self.tpcs.values():
            for ring in ['LER', 'HER']:
                tmp = data[tpc].loc[(data[tpc][ring + '_flag']==1)]# & (data[tpc]['npoints'].isna()==False) & (data[tpc]['Lumi_err'].isna()==False)]
                key = tpc+'_'+ring
                y = array('d', tmp['rate']/tmp['I_'+ring])
                yerr = array('d', tmp['rate']/tmp['I_'+ring]*np.sqrt((tmp['rate_err']/tmp['rate'])**2+(tmp['I_'+ring+'_err']/tmp['I_'+ring])**2))
                x1 = array('d', tmp['P_'+ring])
                x1err = array('d', tmp['P_'+ring+'_err'])
                x2 = array('d', (tmp['I_'+ring]/(tmp['Sy_'+ring]*tmp['Nb_'+ring])))
                x2err = array('d', (tmp['I_'+ring]/(tmp['Sy_'+ring]*tmp['Nb_'+ring]))*np.sqrt((tmp['I_'+ring+'_err']/tmp['I_'+ring])**2+(tmp['Sy_'+ring+'_err']/tmp['Sy_'+ring])**2))
                fits[key + '_B0'], fits[key + '_B1'],fits[key + '_T'],fits[key + '_B0_err'],fits[key + '_B1_err'],fits[key + '_T_err'] = fit(tpc,y,yerr,x1,x1err,x2,x2err)

                #y = array('d', tmp['rate'])
                #yerr = array('d', tmp['rate_err'])
                #x1 = array('d', tmp['I_'+ring]**2)
                #x1err = array('d', (2*tmp['I_'+ring]*tmp['I_'+ring+'_err']))
                #x2 = array('d', (tmp['I_'+ring]**2/(tmp['Sy_'+ring]*tmp['Nb_'+ring])))
                #x2err = array('d', (tmp['I_'+ring]**2/(tmp['Sy_'+ring]*tmp['Nb_'+ring]))*np.sqrt((2*tmp['I_'+ring+'_err']/tmp['I_'+ring])**2+(tmp['Sy_'+ring+'_err']/tmp['Sy_'+ring])**2))
                #fits[key + '_B1'],fits[key + '_T'],fits[key + '_B1_err'],fits[key + '_T_err'] = fit(tpc,y,yerr,x1,x1err,x2,x2err)
                #y = array('d', tmp['rate']/tmp['I_'+ring]**2)
                #yerr = array('d', (tmp['rate']/tmp['I_'+ring]**2)*np.sqrt((tmp['rate_err']/tmp['rate'])**2+(2*tmp['I_'+ring+'_err']/tmp['I_'+ring])**2))
                #x1 = array('d', (1/(tmp['Sy_'+ring]*tmp['Nb_'+ring])))
                #x1err = array('d', (1/(tmp['Sy_'+ring]*tmp['Nb_'+ring]))*(tmp['Sy_'+ring+'_err']/tmp['Sy_'+ring]))
                #fits[key + '_B1'],fits[key + '_T'],fits[key + '_B1_err'],fits[key + '_T_err'] = fit(tpc,y,yerr,x1,x1err)

                try:
                    if fits[key + '_B1_err']/fits[key + '_B1'] > 10:
                        fits[key + '_B1_err']=0
                except ZeroDivisionError:
                    fits[key + '_B1_err']=0
                try:
                    if fits[key + '_B0_err']/fits[key + '_B0'] > 10:
                        fits[key + '_B0_err']=0
                except ZeroDivisionError:
                    fits[key + '_B0_err']=0
                try:
                    if fits[key + '_T_err']/fits[key + '_T'] > 10:
                        fits[key + '_T_err']=0
                except ZeroDivisionError:
                    fits[key + '_T _err']=0
                for fill in ['Cont_inj','Decay']:
                    tmp = data[tpc].loc[(data[tpc][ring + '_flag']==1) & (data[tpc][fill+'_flag']==1)]# & (data[tpc]['npoints'].isna()==False) & (data[tpc]['Lumi_err'].isna()==False)]
                    y = array('d', tmp['rate']/tmp['I_'+ring])
                    yerr = array('d', (tmp['rate']/tmp['I_'+ring])*np.sqrt((tmp['rate_err']/tmp['rate'])**2+(tmp['I_'+ring+'_err']/tmp['I_'+ring])**2))
                    x1 = array('d', tmp['P_'+ring])
                    x1err = array('d', tmp['P_'+ring+'_err'])
                    x2 = array('d', (tmp['I_'+ring]/(tmp['Sy_'+ring]*tmp['Nb_'+ring])))
                    x2err = array('d', (tmp['I_'+ring]/(tmp['Sy_'+ring]*tmp['Nb_'+ring]))*np.sqrt((tmp['I_'+ring+'_err']/tmp['I_'+ring])**2+(tmp['Sy_'+ring+'_err']/tmp['Sy_'+ring])**2))
                    key = tpc+'_'+ring+'_'+fill
                    fits[key + '_B0'], fits[key + '_B1'],fits[key + '_T'],fits[key + '_B0_err'], fits[key + '_B1_err'],fits[key + '_T_err'] = fit(tpc,y,yerr,x1,x1err,x2,x2err)
                    try:
                        if fits[key + '_B1_err']/fits[key + '_B1'] > 10:
                            fits[key + '_B1_err']=0
                    except ZeroDivisionError:
                        fits[key + '_B1_err']=0
                    try:
                        if fits[key + '_B0_err']/fits[key + '_B0'] > 10:
                            fits[key + '_B0_err']=0
                    except ZeroDivisionError:
                        fits[key + '_B0_err']=0
                    try:
                        if fits[key + '_T_err']/fits[key + '_T'] > 10:
                            fits[key + '_T_err']=0
                    except ZeroDivisionError:
                        fits[key + '_T_err']=0
        return fits    

    def apply_fits_to_data(self, period, cont_inj = True, decay = True):
        data = self.binned_data
        fits = self.SB_fit_params
        reduced_data = {}
        for tpc in self.tpcs.values():
            reduced_data[tpc] = data[tpc]
            if (cont_inj and decay) or (cont_inj == False and decay == False):
                reduced_data[tpc] = data[tpc]#.loc[(data[tpc][period+'_flag'] == 1)]
                if period == 'Lumi':
                    reduced_data[tpc]['BG_pred'] = fits[tpc+'_HER_B0']*reduced_data[tpc]['I_HER']*reduced_data[tpc]['P_HER'] + fits[tpc+'_LER_B0']*reduced_data[tpc]['I_LER']*reduced_data[tpc]['P_LER'] - fits[tpc+'_HER_B1']*reduced_data[tpc]['I_HER'] - fits[tpc+'_LER_B1']*reduced_data[tpc]['I_LER'] 
                    reduced_data[tpc]['T_pred'] = 1.2*fits[tpc+'_HER_T']*reduced_data[tpc]['I_HER']**2/(reduced_data[tpc]['Sy_HER']*reduced_data[tpc]['Nb_HER']) + 0.79*fits[tpc+'_LER_T']*reduced_data[tpc]['I_LER']**2/(reduced_data[tpc]['Sy_LER']*reduced_data[tpc]['Nb_LER'])
                    reduced_data[tpc]['BG_pred_err'] = np.sqrt((reduced_data[tpc]['I_HER']*reduced_data[tpc]['P_HER']*fits[tpc+'_HER_B0_err'])**2 + (reduced_data[tpc]['I_HER']*fits[tpc+'_HER_B1_err'])**2 + (fits[tpc+'_HER_B0']*reduced_data[tpc]['I_HER']*reduced_data[tpc]['P_HER_err'])**2 + ((fits[tpc+'_HER_B0']*reduced_data[tpc]['P_HER']-fits[tpc+'_HER_B1'])*reduced_data[tpc]['I_HER_err'])**2 + (reduced_data[tpc]['I_LER']*reduced_data[tpc]['P_LER']*fits[tpc+'_LER_B0_err'])**2 + (reduced_data[tpc]['I_LER']*fits[tpc+'_LER_B1_err'])**2 + (fits[tpc+'_LER_B0']*reduced_data[tpc]['I_LER']*reduced_data[tpc]['P_LER_err'])**2 + ((fits[tpc+'_LER_B0']*reduced_data[tpc]['P_LER']-fits[tpc+'_LER_B1'])*reduced_data[tpc]['I_LER_err'])**2)
                    reduced_data[tpc]['T_pred_err'] = np.sqrt((reduced_data[tpc]['I_HER']**2/(reduced_data[tpc]['Sy_HER']*reduced_data[tpc]['Nb_HER'])*1.2*fits[tpc+'_HER_T_err'])**2 + (2*reduced_data[tpc]['I_HER']*1.2*fits[tpc+'_HER_T']/(reduced_data[tpc]['Sy_HER']*reduced_data[tpc]['Nb_HER'])*reduced_data[tpc]['I_HER_err'])**2 +(reduced_data[tpc]['I_HER']**2*1.2*fits[tpc+'_HER_T']/(reduced_data[tpc]['Sy_HER']**2*reduced_data[tpc]['Nb_HER'])*reduced_data[tpc]['Sy_HER_err'])**2 + (reduced_data[tpc]['I_LER']**2/(reduced_data[tpc]['Sy_LER']*reduced_data[tpc]['Nb_LER'])*0.79*fits[tpc+'_LER_T_err'])**2 + (2*reduced_data[tpc]['I_LER']*0.79*fits[tpc+'_LER_T']/(reduced_data[tpc]['Sy_LER']*reduced_data[tpc]['Nb_LER'])*reduced_data[tpc]['I_LER_err'])**2 +(reduced_data[tpc]['I_LER']**2*0.79*fits[tpc+'_LER_T']/(reduced_data[tpc]['Sy_LER']**2*reduced_data[tpc]['Nb_LER'])*reduced_data[tpc]['Sy_LER_err'])**2)
                else:
                    reduced_data[tpc]['BG_pred'] = fits[tpc+'_'+period+'_B0']*reduced_data[tpc]['I_'+period]*reduced_data[tpc]['P_'+period] - fits[tpc+'_'+period+'_B1']*reduced_data[tpc]['I_'+period]
                    reduced_data[tpc]['T_pred'] = fits[tpc+'_'+period+'_T']*reduced_data[tpc]['I_'+period]**2/(reduced_data[tpc]['Sy_'+period]*reduced_data[tpc]['Nb_'+period])
                    reduced_data[tpc]['BG_pred_err'] = np.sqrt((reduced_data[tpc]['I_'+period]*reduced_data[tpc]['P_'+period]*fits[tpc+'_'+period+'_B0_err'])**2 + (reduced_data[tpc]['I_'+period]*fits[tpc+'_'+period+'_B1_err'])**2 + (fits[tpc+'_'+period+'_B0']*reduced_data[tpc]['I_'+period]*reduced_data[tpc]['P_HER_err'])**2 + ((fits[tpc+'_'+period+'_B0']*reduced_data[tpc]['P_'+period]-fits[tpc+'_'+period+'_B1'])*reduced_data[tpc]['I_HER_err'])**2)
                    reduced_data[tpc]['T_pred_err'] = np.sqrt((reduced_data[tpc]['I_'+period]**2/(reduced_data[tpc]['Sy_'+period]*reduced_data[tpc]['Nb_'+period])*fits[tpc+'_'+period+'_T_err'])**2 + (2*reduced_data[tpc]['I_'+period]*fits[tpc+'_'+period+'_T']/(reduced_data[tpc]['Sy_'+period]*reduced_data[tpc]['Nb_'+period])*reduced_data[tpc]['I_'+period+'_err'])**2 +(reduced_data[tpc]['I_'+period]**2*fits[tpc+'_'+period+'_T']/(reduced_data[tpc]['Sy_'+period]**2*reduced_data[tpc]['Nb_'+period])*reduced_data[tpc]['Sy_'+period+'_err'])**2)
                    
            elif (decay and cont_inj == False):
                reduced_data[tpc] = data[tpc].loc[(data[tpc][period+'_flag'] == 1) & (data[tpc]['Decay_flag']==1)]
                if period == 'Lumi':
                    reduced_data[tpc]['BG_pred'] = fits[tpc+'_HER_Decay_B0']*reduced_data[tpc]['I_HER']*reduced_data[tpc]['P_HER'] + fits[tpc+'_LER_Decay_B0']*reduced_data[tpc]['I_LER']*reduced_data[tpc]['P_LER'] - fits[tpc+'_HER_Decay_B1']*reduced_data[tpc]['I_HER'] - fits[tpc+'_LER_Decay_B1']*reduced_data[tpc]['I_LER']
                    reduced_data[tpc]['T_pred'] = 1.2*fits[tpc+'_HER_Decay_T']*reduced_data[tpc]['I_HER']**2/(reduced_data[tpc]['Sy_HER']*reduced_data[tpc]['Nb_HER']) + 0.79*fits[tpc+'_LER_Decay_T']*reduced_data[tpc]['I_LER']**2/(reduced_data[tpc]['Sy_LER']*reduced_data[tpc]['Nb_LER'])
                    reduced_data[tpc]['BG_pred_err'] = np.sqrt((reduced_data[tpc]['I_HER']*reduced_data[tpc]['P_HER']*fits[tpc+'_HER_Decay_B0_err'])**2 + (reduced_data[tpc]['I_HER']*fits[tpc+'_HER_Decay_B1_err'])**2 + (fits[tpc+'_HER_Decay_B0']*reduced_data[tpc]['I_HER']*reduced_data[tpc]['P_HER_err'])**2 + ((fits[tpc+'_HER_Decay_B0']*reduced_data[tpc]['P_HER']-fits[tpc+'_HER_Decay_B1'])*reduced_data[tpc]['I_HER_err'])**2 + (reduced_data[tpc]['I_LER']*reduced_data[tpc]['P_LER']*fits[tpc+'_LER_Decay_B0_err'])**2 + (reduced_data[tpc]['I_LER']*fits[tpc+'_LER_Decay_B1_err'])**2 + (fits[tpc+'_LER_Decay_B0']*reduced_data[tpc]['I_LER']*reduced_data[tpc]['P_LER_err'])**2 + ((fits[tpc+'_LER_Decay_B0']*reduced_data[tpc]['P_LER']-fits[tpc+'_LER_Decay_B1'])*reduced_data[tpc]['I_LER_err'])**2)
                    reduced_data[tpc]['T_pred_err'] = np.sqrt((reduced_data[tpc]['I_HER']**2/(reduced_data[tpc]['Sy_HER']*reduced_data[tpc]['Nb_HER'])*1.2*fits[tpc+'_HER_Decay_T_err'])**2 + (2*reduced_data[tpc]['I_HER']*1.2*fits[tpc+'_HER_Decay_T']/(reduced_data[tpc]['Sy_HER']*reduced_data[tpc]['Nb_HER'])*reduced_data[tpc]['I_HER_err'])**2 +(reduced_data[tpc]['I_HER']**2*1.2*fits[tpc+'_HER_Decay_T']/(reduced_data[tpc]['Sy_HER']**2*reduced_data[tpc]['Nb_HER'])*reduced_data[tpc]['Sy_HER_err'])**2 + (reduced_data[tpc]['I_LER']**2/(reduced_data[tpc]['Sy_LER']*reduced_data[tpc]['Nb_LER'])*0.79*fits[tpc+'_LER_Decay_T_err'])**2 + (2*reduced_data[tpc]['I_LER']*0.79*fits[tpc+'_LER_Decay_T']/(reduced_data[tpc]['Sy_LER']*reduced_data[tpc]['Nb_LER'])*reduced_data[tpc]['I_LER_err'])**2 +(reduced_data[tpc]['I_LER']**2*0.79*fits[tpc+'_LER_Decay_T']/(reduced_data[tpc]['Sy_LER']**2*reduced_data[tpc]['Nb_LER'])*reduced_data[tpc]['Sy_LER_err'])**2)
                else:
                    reduced_data[tpc]['BG_pred'] = fits[tpc+'_'+period+'_Decay_B0']*reduced_data[tpc]['I_'+period]*reduced_data[tpc]['P_'+period] - fits[tpc+'_'+period+'_Decay_B1']*reduced_data[tpc]['I_'+period]
                    reduced_data[tpc]['T_pred'] = fits[tpc+'_'+period+'_Decay_T']*reduced_data[tpc]['I_'+period]**2/(reduced_data[tpc]['Sy_'+period]*reduced_data[tpc]['Nb_'+period])
                    reduced_data[tpc]['BG_pred_err'] = np.sqrt((reduced_data[tpc]['I_'+period]*reduced_data[tpc]['P_'+period]*fits[tpc+'_'+period+'_Decay_B0_err'])**2 + (reduced_data[tpc]['I_'+period]*fits[tpc+'_'+period+'_Decay_B1_err'])**2 + (fits[tpc+'_'+period+'_Decay_B0']*reduced_data[tpc]['I_'+period]*reduced_data[tpc]['P_HER_err'])**2 + ((fits[tpc+'_'+period+'_Decay_B0']*reduced_data[tpc]['P_'+period]-fits[tpc+'_'+period+'_Decay_B1'])*reduced_data[tpc]['I_HER_err'])**2)
                    reduced_data[tpc]['T_pred_err'] = np.sqrt((reduced_data[tpc]['I_'+period]**2/(reduced_data[tpc]['Sy_'+period]*reduced_data[tpc]['Nb_'+period])*fits[tpc+'_'+period+'_Decay_T_err'])**2 + (2*reduced_data[tpc]['I_'+period]*fits[tpc+'_'+period+'_Decay_T']/(reduced_data[tpc]['Sy_'+period]*reduced_data[tpc]['Nb_'+period])*reduced_data[tpc]['I_'+period+'_err'])**2 +(reduced_data[tpc]['I_'+period]**2*fits[tpc+'_'+period+'_Decay_T']/(reduced_data[tpc]['Sy_'+period]**2*reduced_data[tpc]['Nb_'+period])*reduced_data[tpc]['Sy_'+period+'_err'])**2)
            
            elif (cont_inj and decay == False):
                reduced_data[tpc] = data[tpc].loc[(data[tpc][period+'_flag'] == 1) & (data[tpc]['Cont_inj_flag']==1)]
                if period == 'Lumi':
                    reduced_data[tpc]['BG_pred'] = fits[tpc+'_HER_Decay_B0']*reduced_data[tpc]['I_HER']*reduced_data[tpc]['P_HER'] + fits[tpc+'_LER_Decay_B0']*reduced_data[tpc]['I_LER']*reduced_data[tpc]['P_LER'] - fits[tpc+'_HER_Decay_B1']*reduced_data[tpc]['I_HER'] - fits[tpc+'_LER_Decay_B1']*reduced_data[tpc]['I_LER']
                    reduced_data[tpc]['T_pred'] = 1.2*fits[tpc+'_HER_Decay_T']*reduced_data[tpc]['I_HER']**2/(reduced_data[tpc]['Sy_HER']*reduced_data[tpc]['Nb_HER']) + 0.79*fits[tpc+'_LER_Decay_T']*reduced_data[tpc]['I_LER']**2/(reduced_data[tpc]['Sy_LER']*reduced_data[tpc]['Nb_LER'])
                    reduced_data[tpc]['BG_pred_err'] = np.sqrt((reduced_data[tpc]['I_HER']*reduced_data[tpc]['P_HER']*fits[tpc+'_HER_Decay_B0_err'])**2 + (reduced_data[tpc]['I_HER']*fits[tpc+'_HER_Decay_B1_err'])**2 + (fits[tpc+'_HER_Decay_B0']*reduced_data[tpc]['I_HER']*reduced_data[tpc]['P_HER_err'])**2 + ((fits[tpc+'_HER_Decay_B0']*reduced_data[tpc]['P_HER']-fits[tpc+'_HER_Decay_B1'])*reduced_data[tpc]['I_HER_err'])**2 + (reduced_data[tpc]['I_LER']*reduced_data[tpc]['P_LER']*fits[tpc+'_LER_Decay_B0_err'])**2 + (reduced_data[tpc]['I_LER']*fits[tpc+'_LER_Decay_B1_err'])**2 + (fits[tpc+'_LER_Decay_B0']*reduced_data[tpc]['I_LER']*reduced_data[tpc]['P_LER_err'])**2 + ((fits[tpc+'_LER_Decay_B0']*reduced_data[tpc]['P_LER']-fits[tpc+'_LER_Decay_B1'])*reduced_data[tpc]['I_LER_err'])**2)
                    reduced_data[tpc]['T_pred_err'] = np.sqrt((reduced_data[tpc]['I_HER']**2/(reduced_data[tpc]['Sy_HER']*reduced_data[tpc]['Nb_HER'])*1.2*fits[tpc+'_HER_Decay_T_err'])**2 + (2*reduced_data[tpc]['I_HER']*1.2*fits[tpc+'_HER_Decay_T']/(reduced_data[tpc]['Sy_HER']*reduced_data[tpc]['Nb_HER'])*reduced_data[tpc]['I_HER_err'])**2 +(reduced_data[tpc]['I_HER']**2*1.2*fits[tpc+'_HER_Decay_T']/(reduced_data[tpc]['Sy_HER']**2*reduced_data[tpc]['Nb_HER'])*reduced_data[tpc]['Sy_HER_err'])**2 + (reduced_data[tpc]['I_LER']**2/(reduced_data[tpc]['Sy_LER']*reduced_data[tpc]['Nb_LER'])*0.79*fits[tpc+'_LER_Decay_T_err'])**2 + (2*reduced_data[tpc]['I_LER']*0.79*fits[tpc+'_LER_Decay_T']/(reduced_data[tpc]['Sy_LER']*reduced_data[tpc]['Nb_LER'])*reduced_data[tpc]['I_LER_err'])**2 +(reduced_data[tpc]['I_LER']**2*0.79*fits[tpc+'_LER_Decay_T']/(reduced_data[tpc]['Sy_LER']**2*reduced_data[tpc]['Nb_LER'])*reduced_data[tpc]['Sy_LER_err'])**2)
                else:
                    reduced_data[tpc]['BG_pred'] = fits[tpc+'_'+period+'_Decay_B0']*reduced_data[tpc]['I_'+period]*reduced_data[tpc]['P_'+period] - fits[tpc+'_'+period+'_Decay_B1']*reduced_data[tpc]['I_'+period]
                    reduced_data[tpc]['T_pred'] = fits[tpc+'_'+period+'_Decay_T']*reduced_data[tpc]['I_'+period]**2/(reduced_data[tpc]['Sy_'+period]*reduced_data[tpc]['Nb_'+period])
                    reduced_data[tpc]['BG_pred_err'] = np.sqrt((reduced_data[tpc]['I_'+period]*reduced_data[tpc]['P_'+period]*fits[tpc+'_'+period+'_Decay_B0_err'])**2 + (reduced_data[tpc]['I_'+period]*fits[tpc+'_'+period+'_Decay_B1_err'])**2 + (fits[tpc+'_'+period+'_Decay_B0']*reduced_data[tpc]['I_'+period]*reduced_data[tpc]['P_HER_err'])**2 + ((fits[tpc+'_'+period+'_Decay_B0']*reduced_data[tpc]['P_'+period]-fits[tpc+'_'+period+'_Decay_B1'])*reduced_data[tpc]['I_HER_err'])**2)
                    reduced_data[tpc]['T_pred_err'] = np.sqrt((reduced_data[tpc]['I_'+period]**2/(reduced_data[tpc]['Sy_'+period]*reduced_data[tpc]['Nb_'+period])*fits[tpc+'_'+period+'_Decay_T_err'])**2 + (2*reduced_data[tpc]['I_'+period]*fits[tpc+'_'+period+'_Decay_T']/(reduced_data[tpc]['Sy_'+period]*reduced_data[tpc]['Nb_'+period])*reduced_data[tpc]['I_'+period+'_err'])**2 +(reduced_data[tpc]['I_'+period]**2*fits[tpc+'_'+period+'_Decay_T']/(reduced_data[tpc]['Sy_'+period]**2*reduced_data[tpc]['Nb_'+period])*reduced_data[tpc]['Sy_'+period+'_err'])**2)
            #for i, val in enumerate(reduced_data[tpc]['BG_pred']):
            #    if val <0:
            #        reduced_data[tpc]['T_pred'][i] = reduced_data[tpc]['T_pred'][i] + reduced_data[tpc]['BG_pred'][i]
            #        reduced_data[tpc]['BG_pred'][i] = 0
            
        return reduced_data

    def plot_fits_vs_time(self,ring,tpc,cont_inj=True,decay=True):
        SB_data = self.apply_fits_to_data(ring,cont_inj,decay)
        if (cont_inj and decay) or (cont_inj == False and decay == False):
            key = ring + ' ' + tpc + ' fits to cont_inj and decay'
            lkey = ''
        elif cont_inj == False and decay:
            key = ring + ' ' + tpc + ' fits to decay only'
            lkey = '_Decay'
        elif cont_inj and decay == False:
            key = ring + ' ' + tpc + ' fits to cont_inj only'
            lkey = '_Cont_inj'
        data = {}
        data[tpc] = SB_data[tpc].loc[SB_data[tpc][ring+'_flag']==1]
        plt.errorbar(data[tpc]['ts'],data[tpc]['rate'],data[tpc]['rate_err'],data[tpc]['ts_err'],'o',color='k',label='Data')
        plt.errorbar(data[tpc]['ts'],data[tpc]['BG_pred'],data[tpc]['BG_pred_err'],data[tpc]['ts_err'],'o',color='dodgerblue',label='BG_pred')
        plt.errorbar(data[tpc]['ts'],data[tpc]['T_pred'],data[tpc]['T_pred_err'],data[tpc]['ts_err'],'o',color='magenta',label='T_pred')
        if ring == 'Lumi':
            lumi_fits = self.lumi_fit_params
            data[tpc]['Lumi_pred'] = lumi_fits[tpc+lkey+'_b']+ lumi_fits[tpc+lkey+'_m']*data[tpc]['Lumi']/1e4
            data[tpc]['Lumi_pred_err'] = np.sqrt(lumi_fits[tpc+lkey+'_b_err']**2 + (data[tpc]['Lumi']/1e4*lumi_fits[tpc+lkey+'_m_err'])**2 + (lumi_fits[tpc+lkey+'_m']*data[tpc]['Lumi_err']/1e4)**2)
            plt.errorbar(data[tpc]['ts'],data[tpc]['Lumi_pred'],data[tpc]['Lumi_pred_err'],data[tpc]['ts_err'],'o',color='gold',label='Lumi_pred')
            plt.errorbar(data[tpc]['ts'],data[tpc]['T_pred']+data[tpc]['BG_pred']+data[tpc]['Lumi_pred'],np.sqrt(data[tpc]['T_pred_err']**2 + data[tpc]['BG_pred_err']**2+data[tpc]['Lumi_pred_err']**2),data[tpc]['ts_err'],'o',color='green',label='Total pred')
        else:
            plt.errorbar(data[tpc]['ts'],data[tpc]['T_pred']+data[tpc]['BG_pred'],np.sqrt(data[tpc]['T_pred_err']**2 + data[tpc]['BG_pred_err']**2),data[tpc]['ts_err'],'o',color='green',label='Total pred')
        plt.xticks([])
        plt.xlabel('time')
        plt.ylabel('Rate [Hz]')
        plt.legend(ncol=2)
        plt.title(key)
        plt.tight_layout()
        plt.savefig(ring + '_' + tpc + lkey + '.jpg',dpi=200)
        plt.show()

    def fit_lumi(self):
        l_cont_inj = self.apply_fits_to_data('Lumi', cont_inj = True, decay = False) #can change these to fit assuming decay fit parameters on injection data
        l_decay = self.apply_fits_to_data('Lumi', cont_inj = False, decay = True)
        l = self.apply_fits_to_data('Lumi', cont_inj = True, decay = True)
        lumi_cont_inj = {}
        lumi_decay = {}
        lumi = {}
        lumi_fits = {}

        def fit(tpc,y,yerr,x1,x1err):
            f = ROOT.TF1("f","[0] + [1]*x")
            #f.SetParLimits(0,0,10)
            #f.SetParLimits(1,0,10)
            gr = ROOT.TGraphErrors(len(x1), x1,  y, x1err, yerr)
            gr.Fit(f, 'SEM')
            b = f.GetParameter(0)
            m = f.GetParameter(1)
            b_err = f.GetParError(0)
            m_err = f.GetParError(1)
            return b, m, b_err, m_err

        for tpc in self.tpcs.values():
            lumi_cont_inj[tpc] = l_cont_inj[tpc].loc[(l_cont_inj[tpc]['Lumi_flag']==1) & (l_cont_inj[tpc]['Cont_inj_flag']==1)]
            lumi_decay[tpc] = l_decay[tpc].loc[(l_decay[tpc]['Lumi_flag']==1) & (l_decay[tpc]['Decay_flag']==1)]
            lumi[tpc]= l[tpc].loc[(l[tpc]['Lumi_flag']==1)]
            lumi_cont_inj[tpc]['Scaled_rate'] = lumi_cont_inj[tpc]['rate']-(lumi_cont_inj[tpc]['BG_pred']+lumi_cont_inj[tpc]['T_pred'])
            lumi_cont_inj[tpc]['Scaled_rate_err'] = np.sqrt(lumi_cont_inj[tpc]['rate_err']**2+lumi_cont_inj[tpc]['BG_pred_err']**2+lumi_cont_inj[tpc]['T_pred_err']**2)
            lumi_decay[tpc]['Scaled_rate'] = lumi_decay[tpc]['rate']-(lumi_decay[tpc]['BG_pred']+lumi_decay[tpc]['T_pred'])
            lumi_decay[tpc]['Scaled_rate_err'] = np.sqrt(lumi_decay[tpc]['rate_err']**2+lumi_decay[tpc]['BG_pred_err']**2+lumi_decay[tpc]['T_pred_err']**2)
            lumi[tpc]['Scaled_rate'] = lumi[tpc]['rate']-(lumi[tpc]['BG_pred']+lumi[tpc]['T_pred'])
            lumi[tpc]['Scaled_rate_err'] = np.sqrt(lumi[tpc]['rate_err']**2+lumi[tpc]['BG_pred_err']**2+lumi[tpc]['T_pred_err']**2)
            
            lumi_fits[tpc + '_b'],lumi_fits[tpc + '_m'],lumi_fits[tpc + '_b_err'],lumi_fits[tpc + '_m_err'] = fit(tpc,y=array('d',lumi[tpc]['Scaled_rate']),yerr=array('d',lumi[tpc]['Scaled_rate_err']),x1=array('d',lumi[tpc]['Lumi']/1e4),x1err=array('d',lumi[tpc]['Lumi_err']/1e4))

            lumi_fits[tpc + '_Decay_b'],lumi_fits[tpc + '_Decay_m'],lumi_fits[tpc + '_Decay_b_err'],lumi_fits[tpc + '_Decay_m_err'] = fit(tpc,y=array('d',lumi_decay[tpc]['Scaled_rate']),yerr=array('d',lumi_decay[tpc]['Scaled_rate_err']),x1=array('d',lumi_decay[tpc]['Lumi']/1e4),x1err=array('d',lumi_decay[tpc]['Lumi_err']/1e4))

            lumi_fits[tpc + '_Cont_inj_b'],lumi_fits[tpc + '_Cont_inj_m'],lumi_fits[tpc + '_Cont_inj_b_err'],lumi_fits[tpc + '_Cont_inj_m_err'] = fit(tpc,y=array('d',lumi_cont_inj[tpc]['Scaled_rate']),yerr=array('d',lumi_cont_inj[tpc]['Scaled_rate_err']),x1=array('d',lumi_cont_inj[tpc]['Lumi']/1e4),x1err=array('d',lumi_cont_inj[tpc]['Lumi_err']/1e4))

        return lumi_fits, lumi_cont_inj, lumi_decay, lumi

    def plot_lumi_fits(self):
        fits, cinj, decay, both = self.fit_lumi()
        x = np.linspace(0,2,101)
        plt.figure(figsize=(18,8))
        i=1
        for loc,tpc in zip(self.tpcs.keys(),self.tpcs.values()):
            plt.subplot(2,3,i)
            plt.errorbar(cinj[tpc]['Lumi']/1e4,cinj[tpc]['Scaled_rate'],cinj[tpc]['Scaled_rate_err'],cinj[tpc]['Lumi_err']/1e4,'o',label='Injection',color = 'tab:blue')
            plt.errorbar(decay[tpc]['Lumi']/1e4,decay[tpc]['Scaled_rate'],decay[tpc]['Scaled_rate_err'],decay[tpc]['Lumi_err']/1e4,'o',label='Decay', color = 'tab:orange')
            plt.plot(x,fits[tpc+'_Decay_b']+fits[tpc+'_Decay_m']*x,color='tab:orange')
            plt.plot(x,fits[tpc+'_Cont_inj_b']+fits[tpc+'_Cont_inj_m']*x,color='tab:blue')
            #plt.plot(x,fits[tpc+'_b']+fits[tpc+'_m']*x,color='k',label='Combined fit')
            plt.xlim(0,2)
            #if i <=3:
            #    plt.ylim(-0.1,2.2)
            #else:
            #    plt.ylim(-0.1,0.8)
            plt.title(loc)
            plt.legend()
            plt.xlabel(r'$L$ [$10^{34}\mathrm{cm}^{-2}\mathrm{s}^{-1}$]')
            plt.ylabel(r'$R_L$ [Hz]')
            i+=1
        plt.tight_layout()
        #plt.savefig('Lumi_newest.jpg',dpi=200)
        plt.show()

    def plot_bg_breakdown(self, cont_inj = True, decay = True):
        SB_fit = self.SB_fit_params
        lumi_fit = self.lumi_fit_params
        machine_params = self.data_machine_params
        df = pd.DataFrame()
        if (cont_inj and decay) or (cont_inj == False and decay == False):
            key = ''
            title = 'Continuous injection and decay'
        elif cont_inj and decay == False:
            key = 'Cont_inj_'
            title = 'Continuous injection'
        elif cont_inj == False and decay:
            key = 'Decay_'
            title = 'Decay'
        for tpc in self.tpcs.values():
            LER_bg = SB_fit[tpc+'_LER_'+key+'B0']*machine_params['I_LER']*machine_params['P_LER']-SB_fit[tpc+'_LER_'+key+'B1']*machine_params['I_LER']
            HER_bg = SB_fit[tpc+'_HER_'+key+'B0']*machine_params['I_HER']*machine_params['P_HER']-SB_fit[tpc+'_HER_'+key+'B1']*machine_params['I_HER']
            LER_T = SB_fit[tpc+'_LER_'+key+'T']*machine_params['I_LER']**2/(machine_params['sy_LER']*machine_params['nb_LER'])
            HER_T = SB_fit[tpc+'_HER_'+key+'T']*machine_params['I_HER']**2/(machine_params['sy_HER']*machine_params['nb_HER'])
            lumi = lumi_fit[tpc+'_'+key+'b']+lumi_fit[tpc+'_'+key+'m']*machine_params['lumi']
            total = LER_bg + HER_bg + LER_T + HER_T + lumi
            LER_bg_frac = LER_bg/total*100
            HER_bg_frac = HER_bg/total*100
            LER_T_frac = LER_T/total*100
            HER_T_frac = HER_T/total*100
            lumi_frac = lumi/total*100
            df[tpc] = [LER_bg_frac,LER_T_frac,HER_bg_frac,HER_T_frac,lumi_frac]
        df = df.T
        df.columns = ['LER Beam Gas', 'LER Touschek', 'HER Beam Gas', 'HER Touschek', 'Luminosity']
        colors = ['cyan', 'magenta', 'dodgerblue', 'purple', 'limegreen']
        df.index = self.tpcs.keys()
        fig, ax = plt.subplots(figsize=(16,9))
        df.plot(kind='bar', stacked=True, color = colors, legend = False,ax=ax)
        ax.set_ylabel('Background Fraction [%]')
        ax.set_ylim(0,120)
        plt.xticks(rotation=45, ha="right")
        ax.set_xlabel('TPC z position')
        ax.set_title(title)
        plt.legend(framealpha = 1, ncol = 3)
        plt.tight_layout()
        plt.savefig(key+'summary.jpg',dpi=200)
        plt.show()
        
            

    def get_MC_data(self):
        self.MC_bgTypes = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all', 'RBB_Lumi', 'twoPhoton_Lumi']
        MC = {}
        tree = 'tree_fe4_after_threshold'
        for tpc in self.tpcs.values():
            dir = '~/data/phase3/spring_2020/05-09-20/geant4_simulation/all_events/%s/'%(tpc)
            MC[tpc] = rp.read_root(dir + '%s_combined2.root'%(tpc), key=tree)
            MC[tpc]['bgType'] = MC[tpc]['bgType'].apply(lambda x: self.MC_bgTypes[int(x)])
            if isinstance(self.E_cut,dict):
                MC[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>(self.E_cut[tpc])]
            else:
                MC[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>(self.E_cut)]
            if self.fei4_restrict:
                MC[tpc] = MC[tpc].loc[MC[tpc]['VRC_id'] == 5]
                MC[tpc].index = [i for i in range(0,len(MC[tpc]))]
        return MC


    def get_MC_rates(self): #Scale to luminosity of interest. Units: 1e34cm-2s-1
        tree = 'tree_fe4_after_threshold'
        MC = self.MC_data
        MC_new = {}
        rates = {}
        rates_err = {}
        df = pd.DataFrame()
        df_err = pd.DataFrame()
        lumi_frac = 25/self.MC_machine_params['lumi']
        
        for tpc in self.tpcs.values():
            rates[tpc] = {}
            rates_err[tpc] = {}
            if isinstance(self.E_cut,dict):
                MC_new[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>(self.E_cut[tpc])]
            else:
                MC_new[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>(self.E_cut)]
            if self.recoils_only and self.fei4_restrict:
                MC_new[tpc] = MC_new[tpc].loc[(MC_new[tpc]['is_recoil'] == 1)&(MC_new[tpc]['VRC_id']==5)]
            elif self.recoils_only == False and self.fei4_restrict:
                MC_new[tpc] = MC_new[tpc].loc[MC_new[tpc]['VRC_id']==5]
            elif self.recoils_only and self.fei4_restrict == False:
                MC_new[tpc] = MC_new[tpc].loc[MC_new[tpc]['is_recoil'] == 1]
            else:
                MC_new[tpc] = MC[tpc]
            
            MC_new[tpc].index = [i for i in range(0,len(MC_new[tpc]))]
            for bg in self.MC_bgTypes:
                if tpc == 'iiwi':
                    scale = .9
                elif tpc == 'nene':
                    scale = 0.99005
                elif tpc =='humu':
                    scale = 0.99402
                elif tpc == 'tako':
                    scale = 0.982
                else:
                    scale = 1
                sim_times = {'Brems_HER_dynamic':400.,'Brems_HER_base':400.,'Coulomb_HER_base':40.,
                     'Coulomb_HER_dynamic':40.,'Brems_LER_base':40.,'Brems_LER_dynamic':40.,
                     'Coulomb_LER_base':4.0,'Coulomb_LER_dynamic':4.0,'Touschek_HER_all':1.6,
                     'Touschek_LER_all':0.4,'RBB_Lumi':scale*9.7e-3*lumi_frac,'twoPhoton_Lumi':1e-2*lumi_frac}

                try:
                    rates[tpc][bg] = len(MC_new[tpc].loc[MC_new[tpc]['bgType']==bg])/(sim_times[bg]*100) #100 is from dialing up cross section
                    rates_err[tpc][bg] = np.sqrt(len(MC_new[tpc].loc[MC_new[tpc]['bgType']==bg]))/(sim_times[bg]*100)
                except FileNotFoundError:
                    rates[tpc][bg] = 0
                    rates_err[tpc][bg] = 0
            df = df.append(pd.DataFrame.from_dict(rates[tpc], 'index').T)
            df_err = df_err.append(pd.DataFrame.from_dict(rates_err[tpc], 'index').T)
        df.index = self.tpcs.values()
        df['LER_bg_base'] = (df['Brems_LER_base'] + df['Coulomb_LER_base'])/1200*self.MC_machine_params['I_LER'] #scale by input I/I_ref
        df['HER_bg_base'] = (df['Brems_HER_base'] + df['Coulomb_HER_base'])/1000*self.MC_machine_params['I_HER']
        df['LER_bg_dynamic'] = (df['Brems_LER_dynamic'] + df['Coulomb_LER_dynamic'])/1200**2*self.MC_machine_params['I_LER']**2
        df['HER_bg_dynamic'] = (df['Brems_HER_dynamic'] + df['Coulomb_HER_dynamic'])/1000**2*self.MC_machine_params['I_HER']**2
        df['LER_T'] = df['Touschek_LER_all']/(1200**2/(37*1576))*(self.MC_machine_params['I_LER']**2/(self.MC_machine_params['sy_LER']*self.MC_machine_params['nb_LER']))
        df['HER_T'] = df['Touschek_HER_all']/(1000**2/(36*1576))*(self.MC_machine_params['I_HER']**2/(self.MC_machine_params['sy_HER']*self.MC_machine_params['nb_HER']))
        df['Lumi'] = df['RBB_Lumi'] + df['twoPhoton_Lumi']
        df = df[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
        df_err.index = self.tpcs.values()
        df_err['LER_bg_base'] = (df_err['Brems_LER_base'] + df_err['Coulomb_LER_base'])/1200*self.MC_machine_params['I_LER'] #scale by input I/I_ref
        df_err['HER_bg_base'] = (df_err['Brems_HER_base'] + df_err['Coulomb_HER_base'])/1000*self.MC_machine_params['I_HER']
        df_err['LER_bg_dynamic'] = (df_err['Brems_LER_dynamic'] + df_err['Coulomb_LER_dynamic'])/1200**2*self.MC_machine_params['I_LER']**2
        df_err['HER_bg_dynamic'] = (df_err['Brems_HER_dynamic'] + df_err['Coulomb_HER_dynamic'])/1000**2*self.MC_machine_params['I_HER']**2
        df_err['LER_T'] = df_err['Touschek_LER_all']/(1200**2/(37*1576))*(self.MC_machine_params['I_LER']**2/(self.MC_machine_params['sy_LER']*self.MC_machine_params['nb_LER']))
        df_err['HER_T'] = df_err['Touschek_HER_all']/(1000**2/(36*1576))*(self.MC_machine_params['I_HER']**2/(self.MC_machine_params['sy_HER']*self.MC_machine_params['nb_HER']))
        df_err['Lumi'] = np.sqrt(df_err['RBB_Lumi']**2 + df_err['twoPhoton_Lumi']**2)
        df_err = df_err[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
        df_err.columns = ['LER_bg_base_err', 'LER_bg_dynamic_err', 'LER_T_err', 'HER_bg_base_err', 'HER_bg_dynamic_err', 'HER_T_err', 'Lumi_err']
        df_combined = pd.concat([df,df_err], axis=1)
        return df_combined    

a = analysis(bin_width=60,fit_package='root')
