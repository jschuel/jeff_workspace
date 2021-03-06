import pandas as pd
import uproot as ur
import root_pandas as rp
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import ROOT
from ROOT import TVector3
import array
import matplotlib
from matplotlib.lines import Line2D
from matplotlib import rc
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

rc('text', usetex=False)
pd.set_option('mode.chained_assignment', None) #remove copy warning

class analysis:

    def __init__(self, E_cut = {'palila': 9.0, 'iiwi': 8.8, 'tako': 4.6, 'nene': 5.6, 'elepaio': 6.4, 'humu': 6.4}, E_cut_err = 0, input_file= "/home/jeff/data/phase3/spring_2020/05-09-20/combined_ntuples/05-09_whole_study_even_newerest.root", recoils_only = True, fei4_restrict = True): #enter negative value for E_Cut_err to get low systematic
        self.raw_tpc_data = self.get_tpc_data(recoils_only = recoils_only, E_cut = E_cut, E_cut_err = E_cut_err)
        self.study_data = self.get_raw_study_data(E_cut = E_cut, E_cut_err = 0)
        self.MC_data = self.get_MC_data(E_cut = E_cut, fei4_restrict = fei4_restrict)
        self.MC_rates = self.get_MC_rates(E_cut = E_cut, fei4_restrict = fei4_restrict, recoils_only = recoils_only)
        
    def get_tpc_data(self, input_dir = '/home/jeff/data/phase3/spring_2020/05-09-20/tpc_root_files/', recoils_only = False, E_cut = {'palila': 9.0, 'iiwi': 8.8, 'tako': 4.6, 'nene': 5.6, 'elepaio': 6.4, 'humu': 6.4}, E_cut_err = 0):
        data = {}
        tpcs = ['iiwi', 'humu', 'nene', 'tako', 'palila', 'elepaio']
        for tpc in tpcs:
            #data[tpc] = ur.open(input_dir + "%s_all_newester6.root"%(tpc))[ur.open(input_dir + "%s_all_newester6.root"%(tpc)).keys()[0]].pandas.df(flatten=False)
            data[tpc] = rp.read_root(input_dir + "%s_all_newester8.root"%(tpc))
            ### For flexibility between passing in dictionaries and numbers
            if isinstance(E_cut,dict) and isinstance(E_cut_err,dict):
                data[tpc] = data[tpc].loc[data[tpc]['track_energy']>=(E_cut[tpc]+E_cut_err[tpc])]
            elif isinstance(E_cut,dict) and isinstance(E_cut_err,dict)==False:
                data[tpc] = data[tpc].loc[data[tpc]['track_energy']>=(E_cut[tpc]+E_cut_err)]
            elif isinstance(E_cut,dict)==False and isinstance(E_cut_err,dict):
                data[tpc] = data[tpc].loc[data[tpc]['track_energy']>=(E_cut+E_cut_err[tpc])]
            else:
                data[tpc] = data[tpc].loc[data[tpc]['track_energy']>=(E_cut+E_cut_err)]
            ###
            if recoils_only == True:
                data[tpc] = data[tpc].loc[data[tpc]['is_recoil'] == 1]
            data[tpc]['ts'] = data[tpc]['timestamp_start'].astype('int')
        return data
    
    def get_raw_study_data(self, input_file= "/home/jeff/data/phase3/spring_2020/05-09-20/combined_ntuples/05-09_whole_study_even_newerest.root", E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0):
        #study_data = ur.open(input_file)[ur.open(input_file).keys()[0]].pandas.df(flatten=False)
        study_data = rp.read_root(input_file)
        tpc_data = self.raw_tpc_data
        tpcs = tpc_data.keys()
        dfs = {}
        for tpc in tpcs:
            if isinstance(E_cut,dict) and isinstance(E_cut_err,dict):
                dfs[tpc] = tpc_data[tpc].loc[tpc_data[tpc]['track_energy']>=(E_cut[tpc]+E_cut_err[tpc])]
            elif isinstance(E_cut,dict) and isinstance(E_cut_err,dict)==False:
                dfs[tpc] = tpc_data[tpc].loc[tpc_data[tpc]['track_energy']>=(E_cut[tpc]+E_cut_err)]
            elif isinstance(E_cut,dict)==False and isinstance(E_cut_err,dict):
                dfs[tpc] = tpc_data[tpc].loc[tpc_data[tpc]['track_energy']>=(E_cut+E_cut_err[tpc])]
            else:
                dfs[tpc] = tpc_data[tpc].loc[tpc_data[tpc]['track_energy']>=(E_cut+E_cut_err)]
            dfs[tpc]['ts'] = dfs[tpc]['timestamp_start'].astype('int')
            study_data[tpc+'_neutrons'] = 0
            pv = dfs[tpc].pivot_table(index=['ts'],aggfunc='size')
            pv = pv.loc[pv.index.isin(study_data['ts'])]
            counts = pd.DataFrame()
            counts['ts'] = pv.index
            counts['rate'] = pv.array
            index = study_data.loc[study_data['ts'].astype('int').isin(counts['ts'])].index.to_numpy()
            study_data[tpc+'_neutrons'][index] = counts['rate'].to_numpy()
        return study_data

    def get_MC_data(self, E_cut, fei4_restrict):
        tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all', 'RBB_Lumi', 'twoPhoton_Lumi']
        MC = {}
        tree = 'tree_fe4_after_threshold'
        for tpc in tpcs:
            dir = '~/data/phase3/spring_2020/05-09-20/geant4_simulation/all_events/%s/'%(tpc)
            MC[tpc] = rp.read_root(dir + '%s_combined2.root'%(tpc), key=tree)
            MC[tpc]['bgType'] = MC[tpc]['bgType'].apply(lambda x: bgtype[int(x)])
            if isinstance(E_cut,dict):
                MC[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>(E_cut[tpc])]
            else:
                MC[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>(E_cut)]
            MC[tpc].index = [i for i in range(0,len(MC[tpc]))]
            if fei4_restrict == True:
                MC[tpc] = MC[tpc].loc[MC[tpc]['VRC_id'] == 5]
                MC[tpc].index = [i for i in range(0,len(MC[tpc]))]
        return MC
    
    def get_MC_rates(self, E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0, fei4_restrict = True, recoils_only = True, I_HER = 1000, I_LER = 1200, sy_LER=37, sy_HER=36, nb_LER=1576, nb_HER=1576, lumi=25): #Scale to luminosity of interest. Units: 1e34cm-2s-1
        tpcs = ['elepaio', 'tako', 'palila', 'iiwi', 'nene', 'humu']
        bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all', 'RBB_Lumi', 'twoPhoton_Lumi']
        tree = 'tree_fe4_after_threshold'
        MC = self.MC_data
        MC_new = {}
        rates = {}
        rates_err = {}
        df = pd.DataFrame()
        df_err = pd.DataFrame()
        lumi_frac = 25/lumi
        for tpc in tpcs:
            rates[tpc] = {}
            rates_err[tpc] = {}
            if isinstance(E_cut,dict):
                MC_new[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>(E_cut[tpc])]
            else:
                MC_new[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>(E_cut)]
            if recoils_only == True and fei4_restrict == True:
                MC_new[tpc] = MC_new[tpc].loc[(MC_new[tpc]['is_recoil'] == 1)&(MC_new[tpc]['VRC_id']==5)]
            elif recoils_only == False and fei4_restrict == True:
                MC_new[tpc] = MC_new[tpc].loc[MC_new[tpc]['VRC_id']==5]
            elif recoils_only == True and fei4_restrict == False:
                MC_new[tpc] = MC_new[tpc].loc[MC_new[tpc]['is_recoil'] == 1]
            else:
                MC_new[tpc] = MC[tpc]
            #if fei4_restrict == True:
            #    MC_new[tpc] = MC_new[tpc].loc[MC_new[tpc]['VRC_id']==5] #ID for actual fei4 active area
            MC_new[tpc].index = [i for i in range(0,len(MC_new[tpc]))]
            for bg in bgtype:
                if (bg == 'Brems_HER_dynamic'):
                    t = 400.
                elif (bg == 'Brems_HER_base'):
                    t = 400.
                elif (bg == 'Coulomb_HER_base') or (bg == 'Coulomb_HER_dynamic'):
                    t = 40.
                elif (bg == 'Brems_LER_base') or (bg == 'Brems_LER_dynamic'):
                    t = 40.
                elif (bg == 'Coulomb_LER_base') or (bg == 'Coulomb_LER_dynamic'):
                    t = 4.0
                elif (bg == 'Touschek_HER_all'):
                    t = 1.6
                elif (bg == 'Touschek_LER_all'):
                    t = 0.4
                elif (bg == 'RBB_Lumi'):
                    #t = 2.2e-3*lumi_frac
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
                    t = scale*9.7e-3*lumi_frac
                elif (bg == 'twoPhoton_Lumi'):
                    t = 1e-2*lumi_frac
                try:
                    rates[tpc][bg] = len(MC_new[tpc].loc[MC_new[tpc]['bgType']==bg])/(t*100) #100 is from dialing up cross section
                    rates_err[tpc][bg] = np.sqrt(len(MC_new[tpc].loc[MC_new[tpc]['bgType']==bg]))/(t*100)
                except FileNotFoundError:
                    rates[tpc][bg] = 0
                    rates_err[tpc][bg] = 0
            df = df.append(pd.DataFrame.from_dict(rates[tpc], 'index').T)
            df_err = df_err.append(pd.DataFrame.from_dict(rates_err[tpc], 'index').T)
        df.index = tpcs
        df['LER_bg_base'] = (df['Brems_LER_base'] + df['Coulomb_LER_base'])/1200*I_LER #scale by input I/I_ref
        df['HER_bg_base'] = (df['Brems_HER_base'] + df['Coulomb_HER_base'])/1000*I_HER
        df['LER_bg_dynamic'] = (df['Brems_LER_dynamic'] + df['Coulomb_LER_dynamic'])/1200**2*I_LER**2
        df['HER_bg_dynamic'] = (df['Brems_HER_dynamic'] + df['Coulomb_HER_dynamic'])/1000**2*I_HER**2
        df['LER_T'] = df['Touschek_LER_all']/(1200**2/(37*1576))*(I_LER**2/(sy_LER*nb_LER))
        df['HER_T'] = df['Touschek_HER_all']/(1000**2/(36*1576))*(I_HER**2/(sy_HER*nb_HER))
        df['Lumi'] = df['RBB_Lumi'] + df['twoPhoton_Lumi']
        df = df[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
        df_err.index = tpcs
        df_err['LER_bg_base'] = (df_err['Brems_LER_base'] + df_err['Coulomb_LER_base'])/1200*I_LER #scale by input I/I_ref
        df_err['HER_bg_base'] = (df_err['Brems_HER_base'] + df_err['Coulomb_HER_base'])/1000*I_HER
        df_err['LER_bg_dynamic'] = (df_err['Brems_LER_dynamic'] + df_err['Coulomb_LER_dynamic'])/1200**2*I_LER**2
        df_err['HER_bg_dynamic'] = (df_err['Brems_HER_dynamic'] + df_err['Coulomb_HER_dynamic'])/1000**2*I_HER**2
        df_err['LER_T'] = df_err['Touschek_LER_all']/(1200**2/(37*1576))*(I_LER**2/(sy_LER*nb_LER))
        df_err['HER_T'] = df_err['Touschek_HER_all']/(1000**2/(36*1576))*(I_HER**2/(sy_HER*nb_HER))
        df_err['Lumi'] = np.sqrt(df_err['RBB_Lumi']**2 + df_err['twoPhoton_Lumi']**2)
        df_err = df_err[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
        df_err.columns = ['LER_bg_base_err', 'LER_bg_dynamic_err', 'LER_T_err', 'HER_bg_base_err', 'HER_bg_dynamic_err', 'HER_T_err', 'Lumi_err']
        df_combined = pd.concat([df,df_err], axis=1)
        return df_combined

    def select_study(self, study_type, study_period, E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0): #LER, HER, Lumi, Cont_inj, Decay
        #raw_data = self.study_data
        raw_data = self.get_raw_study_data(E_cut = E_cut, E_cut_err = E_cut_err)
        study_data = raw_data.loc[(raw_data['%s_study_flag'%(study_type)]==1) & (raw_data['%s_flag'%(study_period)] == 1)]
        return study_data

    def get_tpc_data_during_study_period(self, study_type, study_period, E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0):
        study_data = self.select_study(study_type, study_period, E_cut=E_cut, E_cut_err=E_cut_err)
        tpc_data = self.raw_tpc_data
        #tpc_data = self.get_tpc_data(E_cut = E_cut)
        tpcs = tpc_data.keys()
        tpc_study_data = {}
        for tpc in tpcs:
            tpc_study_data[tpc] = tpc_data[tpc].loc[tpc_data[tpc]['ts'].isin(study_data['ts'])]
        return tpc_study_data

    def partition_data_into_subsets(self, study_type, study_period, bins = 6, E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0):
        study_data = self.select_study(study_type, study_period, E_cut=E_cut, E_cut_err = E_cut_err)
        study_data = study_data.reset_index(drop=True)
        partition_indices = [study_data.index.to_list()[0]] + study_data.loc[np.abs(study_data['ts'].diff())>10].index.to_list() + [study_data.index.to_list()[len(study_data)-1]]
        data_subsets = {}
        for i in range(0,len(partition_indices)-1):
            data_subsets['fill_%s'%(i)] = [i for i in range(partition_indices[i], partition_indices[i+1])]
        #print(data_subsets)
        dfs = {}
        for key in data_subsets.keys():
            dfs[key] = np.array_split(study_data.iloc[data_subsets[key]], bins)
        return dfs

    def compute_means_and_errs(self, study_type, study_period,bins = 6,E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0):
        partitioned_data = self.partition_data_into_subsets(study_type, study_period, bins = bins, E_cut=E_cut, E_cut_err = E_cut_err)
        means = pd.DataFrame()
        errs = pd.DataFrame()
        for key in partitioned_data.keys():
            for i in range(0,len(partitioned_data[key])):
                means = means.append(partitioned_data[key][i].apply(lambda x: x.mean()), ignore_index = True)
                errs = errs.append(partitioned_data[key][i].apply(lambda x: x.sem()), ignore_index = True)
        errs.columns = [str(col) + '_err' for col in errs.columns]
        for col in errs.columns:
            means[col] = errs[col]
        means = means.drop(columns = ['LER_study_flag_err', 'HER_study_flag_err', 'Lumi_study_flag_err', 'Cont_inj_flag_err', 'Decay_flag_err', 'Nb_HER_err', 'Nb_LER_err'])
        return means

    def get_fit_parameters(self, study_type, study_period, bins = 6, E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0): #Gives parameters B0, B1, and T defined by Rate/I = B0 + B1*I + T*I/(sy*Nb)
        averaged_data = self.compute_means_and_errs(study_type, study_period,bins = bins, E_cut=E_cut, E_cut_err = E_cut_err)
        tpcs = ['iiwi', 'humu', 'nene', 'tako', 'elepaio', 'palila']
        fit = {}
        averaged_data['heuristic_x'] = averaged_data['I_'+study_type]/(averaged_data['Sy_'+study_type]*averaged_data['Nb_'+study_type])
        averaged_data['heuristic_x_err'] = averaged_data['heuristic_x']*np.sqrt((averaged_data['I_'+study_type+'_err']/averaged_data['I_'+study_type])**2+(averaged_data['Sy_'+study_type+'_err']/averaged_data['Sy_'+study_type])**2)
        for tpc in tpcs:
            
            
            averaged_data[tpc+'_heuristic_y'] = averaged_data[tpc+'_neutrons']/(averaged_data['I_'+study_type])
            averaged_data[tpc+'_heuristic_y_err'] = averaged_data[tpc+'_heuristic_y']*np.sqrt((averaged_data[tpc+'_neutrons_err']/averaged_data[tpc+'_neutrons'])**2+(averaged_data['I_'+study_type+'_err']/averaged_data['I_'+study_type])**2)#+(averaged_data['P_'+study_type+'_err']/averaged_data['P_'+study_type])**2)
            X = averaged_data[['I_'+study_type,'heuristic_x']]
            y = averaged_data['%s_heuristic_y'%(tpc)]
            yerr = averaged_data[tpc+'_heuristic_y_err']
            x1 = array.array('d', averaged_data['I_'+study_type])
            x1err = array.array('d', averaged_data['I_'+study_type+'_err'])
            x2 = array.array('d', averaged_data['heuristic_x'])
            x2err = array.array('d', averaged_data['heuristic_x_err'])
            y_root = array.array('d', y)
            y_root_err = array.array('d', yerr)
            #f2 = ROOT.TF2("f2","[0] + [1]*x + [2]*y", X['I_'+study_type].min(), X['I_'+study_type].max(), X['heuristic_x'].min(), X['heuristic_x'].max())
            f2 = ROOT.TF2("f2","[0] + [1]*x + [2]*y", X['I_'+study_type].min(), X['I_'+study_type].max(), X['heuristic_x'].min(), X['heuristic_x'].max())
            #f2.SetParLimits(0,0,0)
            #f2.SetParLimits(1,0,10)
            #f2.SetParLimits(2,0,10)
            #f2.SetParLimits(0,0,1e-5)
            #f2.SetParLimits(1,0,1e-7)
            #f2.SetParLimits(2,0,1e-1)
            
            gr = ROOT.TGraph2DErrors(len(x1), x1, x2, y_root, x1err, x2err, y_root_err)
            gr.Fit(f2, 'SREM')
            fit[tpc+'_B0'] = f2.GetParameter(0)
            fit[tpc+'_B1'] = f2.GetParameter(1)
            fit[tpc+'_T'] = f2.GetParameter(2)
            fit[tpc+'_B0_err'] = f2.GetParError(0)
            fit[tpc+'_B1_err'] = f2.GetParError(1)
            fit[tpc+'_T_err'] = f2.GetParError(2)
            
            '''
            y = array.array('d', averaged_data[tpc+'_neutrons'])
            yerr = array.array('d', averaged_data[tpc+'_neutrons_err'])
            x1 = array.array('d', averaged_data['I_'+study_type])
            x1err = array.array('d', averaged_data['I_'+study_type+'_err'])
            x2 = array.array('d', averaged_data['I_'+study_type]**2/(averaged_data['Sy_'+study_type]*averaged_data['Nb_'+study_type]))
            T_err = np.sqrt((2*averaged_data['I_'+study_type]/(averaged_data['Sy_'+study_type]*averaged_data['Nb_'+study_type])*averaged_data['I_'+study_type+'_err'])**2
                            + (averaged_data['I_'+study_type]**2/(averaged_data['Sy_'+study_type]**2*averaged_data['Nb_'+study_type])*averaged_data['Sy_'+study_type+'_err'])**2)
            x2err = array.array('d', T_err)
            f2 = ROOT.TF2("f2","[0]*x + [1]*x^2 + [2]*y + [3]", 0, averaged_data['I_'+study_type].max(), 0, (averaged_data['I_'+study_type]**2/(averaged_data['Sy_'+study_type]*averaged_data['Nb_'+study_type])).max())
            f2.SetParLimits(0,0,1)
            f2.SetParLimits(1,0,1)
            f2.SetParLimits(2,0,1)
            f2.SetParLimits(3,0,1)
            gr = ROOT.TGraph2DErrors(len(x1), x1, x2, y, x1err, x2err, yerr)
            gr.Fit(f2, 'SM')
            fit[tpc+'_B0'] = f2.GetParameter(0)
            fit[tpc+'_B1'] = f2.GetParameter(1)
            fit[tpc+'_T'] = f2.GetParameter(2)
            fit[tpc+'_D'] = f2.GetParameter(3)
            fit[tpc+'_B0_err'] = f2.GetParError(0)
            fit[tpc+'_B1_err'] = f2.GetParError(1)
            fit[tpc+'_T_err'] = f2.GetParError(2)
            fit[tpc+'_D_err'] = f2.GetParError(3)
            '''
        return fit

    def measure_and_fit_lumi_bgs(self, study_period,bins = 15, E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0):
        if study_period.lower() == "decay":
            bins = 6
        lumi_data_avg = self.compute_means_and_errs("Lumi", study_period,bins = bins, E_cut = E_cut, E_cut_err = E_cut_err)
        LER_fit_params = self.get_fit_parameters("LER", study_period, bins=bins, E_cut = E_cut, E_cut_err = E_cut_err)
        HER_fit_params = self.get_fit_parameters("HER", study_period,bins = bins, E_cut = E_cut, E_cut_err = E_cut_err)
        lumi = self.select_study("Lumi", study_period)
        ler = self.select_study("LER", study_period)
        her = self.select_study("HER", study_period)
        #ler_scale = 1/(lumi['Sy_LER'].mean()/ler['Sy_LER'].mean()) #scale touschek contriubtions down to luminosity conditions (for beam-beam blowup correction)
        ler_scale = 1
        #her_scale = 1/(lumi['Sy_HER'].mean()/her['Sy_HER'].mean())
        her_scale = 1
        tpcs = ['iiwi', 'humu', 'nene', 'tako', 'elepaio', 'palila']
        fits = {}
        fits_corrected = {}
        LER_rates = {}
        HER_rates = {}
        LER_rates_scale = {}
        HER_rates_scale = {}
        LER_rates_err = {}
        HER_rates_err = {}
        lumi_rates = {}
        lumi_rates_scale = {}
        lumi_rates_err = {}
        for tpc in tpcs:
            LER_rates[tpc] = LER_fit_params[tpc+'_B0']*lumi_data_avg['I_LER'] + LER_fit_params[tpc+'_B1']*lumi_data_avg['I_LER']**2 + LER_fit_params[tpc+'_T']*lumi_data_avg['I_LER']**2/(lumi_data_avg['Sy_LER']*lumi_data_avg['Nb_LER'])
            LER_rates_scale[tpc] = LER_fit_params[tpc+'_B0']*lumi_data_avg['I_LER'] + LER_fit_params[tpc+'_B1']*lumi_data_avg['I_LER']**2 + ler_scale**2*LER_fit_params[tpc+'_T']*lumi_data_avg['I_LER']**2/(lumi_data_avg['Sy_LER']*lumi_data_avg['Nb_LER'])
            HER_rates[tpc] = HER_fit_params[tpc+'_B0']*lumi_data_avg['I_HER'] + HER_fit_params[tpc+'_B1']*lumi_data_avg['I_HER']**2 + HER_fit_params[tpc+'_T']*lumi_data_avg['I_HER']**2/(lumi_data_avg['Sy_HER']*lumi_data_avg['Nb_HER'])
            HER_rates_scale[tpc] = HER_fit_params[tpc+'_B0']*lumi_data_avg['I_HER'] + HER_fit_params[tpc+'_B1']*lumi_data_avg['I_HER']**2 + her_scale**2*HER_fit_params[tpc+'_T']*lumi_data_avg['I_HER']**2/(lumi_data_avg['Sy_HER']*lumi_data_avg['Nb_HER'])
            LER_rates_err[tpc] = np.sqrt((LER_fit_params[tpc+'_B0']+2*LER_fit_params[tpc+'_B1']*lumi_data_avg['I_LER']+2*LER_fit_params[tpc+'_T']*lumi_data_avg['I_LER']/(lumi_data_avg['Sy_LER']**2*lumi_data_avg['Nb_LER'])*lumi_data_avg['I_LER_err'])**2+(LER_fit_params[tpc+'_T']*lumi_data_avg['I_LER']**2/(lumi_data_avg['Sy_LER']**2*lumi_data_avg['Nb_LER'])*lumi_data_avg['Sy_LER_err'])**2)
            HER_rates_err[tpc] = np.sqrt((HER_fit_params[tpc+'_B0']+2*HER_fit_params[tpc+'_B1']*lumi_data_avg['I_HER']+2*HER_fit_params[tpc+'_T']*lumi_data_avg['I_HER']/(lumi_data_avg['Sy_HER']**2*lumi_data_avg['Nb_HER'])*lumi_data_avg['I_HER_err'])**2+(HER_fit_params[tpc+'_T']*lumi_data_avg['I_HER']**2/(lumi_data_avg['Sy_HER']**2*lumi_data_avg['Nb_HER'])*lumi_data_avg['Sy_HER_err'])**2)
            lumi_rates[tpc] = lumi_data_avg[tpc+'_neutrons'] - LER_rates[tpc] - HER_rates[tpc]
            lumi_rates_scale[tpc] = lumi_data_avg[tpc+'_neutrons'] - LER_rates_scale[tpc] - HER_rates_scale[tpc]
            lumi_rates_err[tpc] = np.sqrt(lumi_data_avg[tpc+'_neutrons_err']**2 + LER_rates_err[tpc]**2 + HER_rates_err[tpc]**2)
        
            gr = ROOT.TGraphErrors(len(lumi_rates[tpc]), array.array('d', lumi_data_avg['ECL_lumi']/10000), array.array('d', lumi_rates[tpc]), array.array('d', lumi_data_avg['ECL_lumi_err']/10000), array.array('d', lumi_rates_err[tpc]))
            gr_corrected = ROOT.TGraphErrors(len(lumi_rates_scale[tpc]), array.array('d', lumi_data_avg['ECL_lumi']/10000), array.array('d', lumi_rates_scale[tpc]), array.array('d', lumi_data_avg['ECL_lumi_err']/10000), array.array('d', lumi_rates_err[tpc]))
            f1 = ROOT.TF1("f1", "[0] + [1]*x", 0, 2)
            gr.Fit("f1", "SEMR")
            fits['%s_int'%(tpc)] = gr.GetFunction("f1").GetParameter(0)
            fits['%s_int_err'%(tpc)] = gr.GetFunction("f1").GetParError(0)
            fits['%s_slope'%(tpc)] = gr.GetFunction("f1").GetParameter(1)
            fits['%s_slope_err'%(tpc)] = gr.GetFunction("f1").GetParError(1)
            f2 = ROOT.TF1("f2", "[0] + [1]*x", 0, 2)
            gr_corrected.Fit("f2", "SEMR")
            fits_corrected['%s_int'%(tpc)] = gr_corrected.GetFunction("f2").GetParameter(0)
            fits_corrected['%s_int_err'%(tpc)] = gr_corrected.GetFunction("f2").GetParError(0)
            fits_corrected['%s_slope'%(tpc)] = gr_corrected.GetFunction("f2").GetParameter(1)
            fits_corrected['%s_slope_err'%(tpc)] = gr_corrected.GetFunction("f2").GetParError(1)
            
        return fits, lumi_rates, lumi_rates_err, lumi_rates_scale, fits_corrected

    def plot_fit(self, lumi_only=True, tunnel = 'BWD', legend = True):
        plt.rc('legend', fontsize=13)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        if lumi_only == True:
            fig, (ax,ax0) = plt.subplots(2,1, figsize = (18,7), gridspec_kw={'height_ratios': [3,1]}, sharex='all')
        else:
            fig, ax = plt.subplots(figsize = (18,9))
        ax1 = ax.twinx()
        ax.set_ylabel(r'Lumi. [$10^{32}$cm$^{-2}$s$^{-1}$], Current [mA]')
        ax.set_ylim(0,600)
        if lumi_only == True:
            ax.set_xlim(-0.03,2.18)
        ax.set_xlabel('Elapsed Time [h]')

        if legend == True:
            if tunnel.lower() == 'bwd':
                labels = ['\n z=-14m','\n z=-8.0m','\n z=-5.6m']
                labels_bot = ['z=-14m','z=-8.0m','z=-5.6m']
            else:
                labels = ['\n z=+6.6m','\n z=+14m','\n z=+16m']
                labels_bot = ['z=+6.6m','z=+14m','z=+16m']
            #shapes = ['o','s','^']
            #colors_data = ['black', 'dimgray', 'silver']
            #colors = ['darkgreen', 'forestgreen', 'lime']
            shapes = ['o','s','^']
            colors_data = ['indigo', 'darkviolet', 'magenta']
            colors = ['darkgreen', 'forestgreen', 'springgreen']
            skb_handles = [Line2D([0], [0], color='b', lw=4, label='HER Current'),Line2D([0], [0], color='r', lw=4, label='LER Current'), Line2D([0], [0], marker='o', color='w', label='Luminosity',markerfacecolor='gold', markersize=15), Line2D([0], [0], marker=shapes[0], color='w', label=labels[0], markerfacecolor=colors_data[0], markersize=0), Line2D([0], [0], marker=shapes[0], color='w', label='Measured', markerfacecolor=colors_data[0], markersize=15), Line2D([0], [0], marker=shapes[0], color='w', label='Fit', markerfacecolor=colors[0], markersize=15), Line2D([0], [0], marker=shapes[1], color='w', label=labels[1], markerfacecolor=colors_data[1], markersize=0), Line2D([0], [0], marker=shapes[1], color='w', label='Measured', markerfacecolor=colors_data[1], markersize=15), Line2D([0], [0], marker=shapes[1], color='w', label='Fit', markerfacecolor=colors[1], markersize=15), Line2D([0], [0], marker=shapes[2], color='w', label=labels[2], markerfacecolor=colors_data[0], markersize=0), Line2D([0], [0], marker=shapes[2], color='w', label='Measured', markerfacecolor=colors_data[2], markersize=15), Line2D([0], [0], marker=shapes[2], color='w', label='Fit', markerfacecolor=colors[2], markersize=15)]
            l_skb = plt.legend(handles = skb_handles,loc = 'upper left', ncol = 1, bbox_to_anchor = (1.1,1))
            
            #l_tpc = plt.legend(handles = tpc_handles,ncol=3, loc='upper center')
            #ax.add_artist(l_skb)
            #ax.add_artist(l_tpc)
        
        ax1.set_ylabel('Rate [Hz]',rotation = 270,labelpad = 30)
        if lumi_only == False:
            ax1.set_ylim(2e-3,20)
            ax1.set_yscale("Log")
        else:
            if tunnel.upper() == 'BWD':
                ax1.set_ylim(0,3)
                ax1.grid()
            else:
                ax1.set_ylim(0,0.6)
                ax1.grid()
        fit_params = {}
        data = {}
        data_avg = {}
        fit_bg_base = {}
        fit_bg_base_err = {}
        fit_bg_dynamic = {}
        fit_bg_dynamic_err = {}
        fit_t = {}
        fit_t_err = {}
        fit = {}
        fit_err = {}
        fit_bg_base_avg = {}
        fit_bg_base_avg_err = {}
        fit_bg_dynamic_avg = {}
        fit_bg_dynamic_avg_err = {}
        fit_t_avg = {}
        fit_t_avg_err = {}
        fit_avg = {}
        fit_avg_err = {}
        if lumi_only == False:
            t0 = self.compute_means_and_errs("LER", "Cont_inj", bins = 6)['ts'][0]
        else:
            t0 = self.compute_means_and_errs("Lumi", "Cont_inj", bins = 15)['ts'][0]
        for study_period in ["Cont_inj", "Decay"]:
            if study_period == "Decay":
                nbins = 6
            else:
                nbins = 15
            fit_params[study_period+'_Lumi'] = self.measure_and_fit_lumi_bgs(study_period, bins = nbins)[4]
            data[study_period+'_Lumi'] = self.select_study('Lumi', study_period)
            data_avg[study_period+'_Lumi'] = self.compute_means_and_errs('Lumi', study_period, bins = nbins)
            ax.plot((data[study_period+'_'+'Lumi']['ts']-t0)/3600, data[study_period+'_'+'Lumi']['I_LER'], 'o', markersize = 1, color = 'red', label = "I_LER [mA]")
            ax.plot((data[study_period+'_'+'Lumi']['ts']-t0)/3600, data[study_period+'_'+'Lumi']['I_HER'], 'o', markersize = 1, color = 'blue', label = "I_HER [mA]")
            ax.plot((data[study_period+'_'+'Lumi']['ts']-t0)/3600, data[study_period+'_'+'Lumi']['ECL_lumi']/100, 'o', markersize = 1, color = 'gold', label = 'Luminosity [a.u.]')

            
            for ring in ['LER', 'HER']:
                fit_params[study_period+'_'+ring] = self.get_fit_parameters(ring, study_period, bins = 6)
                data[study_period+'_'+ring] = self.select_study(ring, study_period)
                data_avg[study_period+'_'+ring] = self.compute_means_and_errs(ring, study_period, bins = 6)
                if lumi_only == False:
                    if ring == 'LER':
                        ax.plot((data[study_period+'_'+ring]['ts']-t0)/3600, data[study_period+'_'+ring]['I_%s'%(ring)], 'o', markersize = 1, color = 'red', label = "I_LER [mA]")
                    else:
                        ax.plot((data[study_period+'_'+ring]['ts']-t0)/3600, data[study_period+'_'+ring]['I_%s'%(ring)], 'o', markersize = 1, color = 'blue', label = "I_HER [mA]")
            for ring in ['LER','HER']:
                #scale_ler = data[study_period+'_'+'LER']['Sy_LER'].mean()/data[study_period+'_'+'Lumi']['Sy_LER'].mean() #scale touschek contriubtions down to luminosity conditions (for beam-beam blowup correction)
                scale_ler = 1
                #scale_her = data[study_period+'_'+'HER']['Sy_HER'].mean()/data[study_period+'_'+'Lumi']['Sy_HER'].mean() #scale touschek contriubtions down to luminosity conditions (for beam-beam blowup correction)
                scale_her = 1
                scale = 1
                shapes = ['o','s','^']
                colors_data = ['indigo', 'darkviolet', 'magenta']
                colors = ['darkgreen', 'forestgreen', 'springgreen']
                i=0
                offsets = [-0.01,0,0.01]
                if tunnel.lower() == 'bwd':
                    tpcs = ['elepaio', 'tako', 'palila']
                else:
                    tpcs = ['iiwi', 'nene', 'humu']
                for tpc in tpcs:
                    fit_bg_base[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_B0']*data[study_period+'_'+ring]['I_%s'%(ring)]
                    fit_bg_dynamic[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_B1']*data[study_period+'_'+ring]['I_%s'%(ring)]**2
                    fit_bg_base_err[study_period+'_'+ring] = np.sqrt((fit_params[study_period+'_'+ring][tpc+'_B0_err']*data[study_period+'_'+ring]['I_%s'%(ring)])**2)
                    fit_bg_dynamic_err[study_period+'_'+ring] = np.sqrt((fit_params[study_period+'_'+ring][tpc+'_B1_err']*data[study_period+'_'+ring]['I_%s'%(ring)]**2)**2)
                    fit_t[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_T']*data[study_period+'_'+ring]['I_%s'%(ring)]**2/(data[study_period+'_'+ring]['Sy_%s'%(ring)]*data[study_period+'_'+ring]['Nb_%s'%(ring)])
                    fit_t_err[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_T_err']*data[study_period+'_'+ring]['I_%s'%(ring)]**2/(data[study_period+'_'+ring]['Sy_%s'%(ring)]*data[study_period+'_'+ring]['Nb_%s'%(ring)])
                    fit[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_B0']*data[study_period+'_'+ring]['I_%s'%(ring)] + fit_params[study_period+'_'+ring][tpc+'_B1']*data[study_period+'_'+ring]['I_%s'%(ring)]**2 + fit_params[study_period+'_'+ring][tpc+'_T']*data[study_period+'_'+ring]['I_%s'%(ring)]**2/(data[study_period+'_'+ring]['Sy_%s'%(ring)]*data[study_period+'_'+ring]['Nb_%s'%(ring)])
                    fit_err[study_period+'_'+ring] = np.sqrt(fit_bg_base_err[study_period+'_'+ring]**2 + fit_bg_dynamic_err[study_period+'_'+ring]**2 + fit_t_err[study_period+'_'+ring]**2)
            
                    fit_bg_base_avg[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_B0']*data_avg[study_period+'_'+ring]['I_%s'%(ring)]
                    fit_bg_dynamic_avg[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_B1']*data_avg[study_period+'_'+ring]['I_%s'%(ring)]**2
                    fit_bg_base_avg_err[study_period+'_'+ring] = np.sqrt((fit_params[study_period+'_'+ring][tpc+'_B0_err']*data_avg[study_period+'_'+ring]['I_%s'%(ring)])**2)
                    fit_bg_dynamic_avg_err[study_period+'_'+ring] = np.sqrt((fit_params[study_period+'_'+ring][tpc+'_B1_err']*data_avg[study_period+'_'+ring]['I_%s'%(ring)]**2)**2)
                    fit_t_avg[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_T']*data_avg[study_period+'_'+ring]['I_%s'%(ring)]**2/(data_avg[study_period+'_'+ring]['Sy_%s'%(ring)]*data_avg[study_period+'_'+ring]['Nb_%s'%(ring)])
                    fit_t_avg_err[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_T_err']*data_avg[study_period+'_'+ring]['I_%s'%(ring)]**2/(data_avg[study_period+'_'+ring]['Sy_%s'%(ring)]*data_avg[study_period+'_'+ring]['Nb_%s'%(ring)])
                    fit_avg[study_period+'_'+ring] = fit_params[study_period+'_'+ring][tpc+'_B0']*data_avg[study_period+'_'+ring]['I_%s'%(ring)] + fit_params[study_period+'_'+ring][tpc+'_B1']*data_avg[study_period+'_'+ring]['I_%s'%(ring)]**2 + scale**2*fit_params[study_period+'_'+ring][tpc+'_T']*data_avg[study_period+'_'+ring]['I_%s'%(ring)]**2/(data_avg[study_period+'_'+ring]['Sy_%s'%(ring)]*data_avg[study_period+'_'+ring]['Nb_%s'%(ring)])
                    fit_avg_err[study_period+'_'+ring] = np.sqrt(fit_bg_base_avg_err[study_period+'_'+ring]**2 + fit_bg_dynamic_avg_err[study_period+'_'+ring]**2 + fit_t_avg_err[study_period+'_'+ring]**2)
                    if lumi_only == False:
                        p1 = ax1.errorbar((data_avg[study_period+'_'+ring]['ts']-t0)/3600, data_avg[study_period+'_'+ring][tpc+'_neutrons'], data_avg[study_period+'_'+ring][tpc+'_neutrons_err'], data_avg[study_period+'_'+ring]['ts_err']/3600, shapes[i], markersize = 6, color = colors_data[i], markeredgecolor = 'k', label = 'data', elinewidth=0.5, alpha = 0.5)
                        p2 = ax1.errorbar((data_avg[study_period+'_'+ring]['ts']-t0)/3600, fit_avg[study_period+'_'+ring], fit_avg_err[study_period+'_'+ring], data_avg[study_period+'_'+ring]['ts_err']/3600, shapes[i], markersize = 6, color = colors[i], markeredgecolor = 'k', label = 'total fit', elinewidth=0.5,alpha = 0.5)
                        p3 = ax1.errorbar((data_avg[study_period+'_'+ring]['ts']-t0)/3600, fit_bg_base_avg[study_period+'_'+ring], fit_bg_base_avg_err[study_period+'_'+ring], data_avg[study_period+'_'+ring]['ts_err']/3600, shapes[i], markersize = 6, color = 'purple', markeredgecolor = 'k', label = 'beam gas base', alpha = 0.6)
                        p4 = ax1.errorbar((data_avg[study_period+'_'+ring]['ts']-t0)/3600, fit_bg_dynamic_avg[study_period+'_'+ring], fit_bg_dynamic_avg_err[study_period+'_'+ring], data_avg[study_period+'_'+ring]['ts_err']/3600, shapes[i], markersize = 6, color = 'gray', markeredgecolor = 'k', label = 'beam gas dyn.', alpha = 0.6)
                        p5 = ax1.errorbar((data_avg[study_period+'_'+ring]['ts']-t0)/3600, fit_t_avg[study_period+'_'+ring], fit_t_avg_err[study_period+'_'+ring], data_avg[study_period+'_'+ring]['ts_err']/3600, shapes[i], markersize = 6, color = 'green', markeredgecolor = 'k', label = 'Touschek', alpha = 0.6)

                    LER_rates = fit_params[study_period+'_'+'LER'][tpc+'_B0']*data_avg[study_period+'_Lumi']['I_LER'] + fit_params[study_period+'_'+'LER'][tpc+'_B1']*data_avg[study_period+'_Lumi']['I_LER']**2 + scale_ler**2*fit_params[study_period+'_'+'LER'][tpc+'_T']*data_avg[study_period+'_Lumi']['I_LER']**2/(data_avg[study_period+'_Lumi']['Sy_LER']*data_avg[study_period+'_Lumi']['Nb_LER']) #during lumi period
                    HER_rates = fit_params[study_period+'_'+'HER'][tpc+'_B0']*data_avg[study_period+'_Lumi']['I_HER'] + fit_params[study_period+'_'+'HER'][tpc+'_B1']*data_avg[study_period+'_Lumi']['I_HER']**2 + scale_her**2*fit_params[study_period+'_'+'HER'][tpc+'_T']*data_avg[study_period+'_Lumi']['I_HER']**2/(data_avg[study_period+'_Lumi']['Sy_HER']*data_avg[study_period+'_Lumi']['Nb_HER']) #during lumi period
                    LER_rates_err = np.sqrt((fit_params[study_period+'_'+'LER'][tpc+'_B0']+2*fit_params[study_period+'_'+'LER'][tpc+'_B1']*data_avg[study_period+'_Lumi']['I_LER']+2*fit_params[study_period+'_'+'LER'][tpc+'_T']*data_avg[study_period+'_Lumi']['I_LER']/(data_avg[study_period+'_Lumi']['Sy_LER']**2*data_avg[study_period+'_Lumi']['Nb_LER'])*data_avg[study_period+'_Lumi']['I_LER_err'])**2+(fit_params[study_period+'_'+'LER'][tpc+'_T']*data_avg[study_period+'_Lumi']['I_LER']**2/(data_avg[study_period+'_Lumi']['Sy_LER']**2*data_avg[study_period+'_Lumi']['Nb_LER'])*data_avg[study_period+'_Lumi']['Sy_LER_err'])**2) #during lumi period
                    HER_rates_err = np.sqrt((fit_params[study_period+'_'+'HER'][tpc+'_B0']+2*fit_params[study_period+'_'+'HER'][tpc+'_B1']*data_avg[study_period+'_Lumi']['I_HER']+2*fit_params[study_period+'_'+'HER'][tpc+'_T']*data_avg[study_period+'_Lumi']['I_HER']/(data_avg[study_period+'_Lumi']['Sy_HER']**2*data_avg[study_period+'_Lumi']['Nb_HER'])*data_avg[study_period+'_Lumi']['I_HER_err'])**2+(fit_params[study_period+'_'+'HER'][tpc+'_T']*data_avg[study_period+'_Lumi']['I_HER']**2/(data_avg[study_period+'_Lumi']['Sy_HER']**2*data_avg[study_period+'_Lumi']['Nb_HER'])*data_avg[study_period+'_Lumi']['Sy_HER_err'])**2) #during lumi period
                    Lumi_rates = LER_rates+HER_rates+fit_params[study_period+'_'+'Lumi'][tpc+'_int']+fit_params[study_period+'_'+'Lumi'][tpc+'_slope']*data_avg[study_period+'_Lumi']['ECL_lumi']/10000
                    Lumi_fit_err = np.sqrt(fit_params[study_period+'_'+'Lumi'][tpc+'_int_err']**2+(data_avg[study_period+'_Lumi']['ECL_lumi']/10000)**2*fit_params[study_period+'_'+'Lumi'][tpc+'_slope_err']**2+fit_params[study_period+'_'+'Lumi'][tpc+'_slope']**2*(data_avg[study_period+'_Lumi']['ECL_lumi_err']/10000)**2)
                    Lumi_rates_err = np.sqrt(LER_rates_err**2 + HER_rates_err**2 + Lumi_fit_err**2)
            
                    ax1.errorbar((data_avg[study_period+'_'+'Lumi']['ts']-t0)/3600, Lumi_rates,Lumi_rates_err, data_avg[study_period+'_Lumi']['ts_err']/3600, shapes[i], markersize = 8, color = colors[i], markeredgecolor = 'k', label = 'Fit', alpha = 0.5, elinewidth = 1)
                    ax1.errorbar((data_avg[study_period+'_Lumi']['ts']-t0)/3600, data_avg[study_period+'_Lumi'][tpc+'_neutrons'], data_avg[study_period+'_Lumi'][tpc+'_neutrons_err'], data_avg[study_period+'_Lumi']['ts_err']/3600, shapes[i], markersize = 8, color = colors_data[i], markeredgecolor = 'k', label = 'data',alpha = 0.5, elinewidth=1)
                    if lumi_only == True:
                        residual = ((data_avg[study_period+'_Lumi'][tpc+'_neutrons'] - Lumi_rates)/data_avg[study_period+'_Lumi'][tpc+'_neutrons'])*100
                        residual_err = np.sqrt((Lumi_rates_err)**2 + (data_avg[study_period+'_Lumi'][tpc+'_neutrons_err'])**2)/(data_avg[study_period+'_Lumi'][tpc+'_neutrons'])*100
                        ax0.errorbar((data_avg[study_period+'_'+'Lumi']['ts']-t0)/3600+offsets[i], residual, residual_err, data_avg[study_period+'_Lumi']['ts_err']/3600, shapes[i], color = 'k', label = labels[i], markeredgecolor='k', markerfacecolor = 'white', alpha = 0.5, markersize = 8)
                    i+=1
        if lumi_only == True:
            ax0.plot(np.linspace(-0.1,2.5,101), np.linspace(0,0,101), '--', color = 'red')
            ax0.set_ylabel(r'$\frac{\mathrm{Meas.} - \mathrm{Fit}}{\mathrm{Meas.}}$ [%]')
            ax0.set_xlabel('Elapsed Time [h]')
            ax0.set_ylim(-80,80)
            ax0.set_yticks([-40,0,40])
            ax0.set_xlim(-0.03,2.18)
            handles = [Line2D([0], [0], marker=shapes[0], color='w', label=labels_bot[0], markeredgecolor='k', markerfacecolor = None, markersize=12, alpha = 0.5), Line2D([0], [0], marker=shapes[1], color='w', label=labels_bot[1], markeredgecolor='k', markerfacecolor = None, markersize=12, alpha = 0.5), Line2D([0], [0], marker=shapes[2], color='w', label=labels_bot[2], markeredgecolor='k', markerfacecolor = None, markersize=12, alpha = 0.5)]
            if tunnel.upper() == 'BWD':
                ax0.legend(loc = 'upper left', handles=handles, ncol=3)
                ax0.text(1.72,26,'(ii)',size = 24)
                ax.text(1.72,450,'(i)',size=24)
            else:
                ax0.legend(loc = 'lower left', handles=handles, ncol=3)
                ax0.text(1.72,26,'(iv)',size = 24)
                ax.text(1.72,450,'(iii)',size=24)
            ax0.yaxis.set_minor_locator(AutoMinorLocator())
            ax0.grid(which = 'both', axis = 'y')
            #ax0.set_title('Residuals')
        plt.subplots_adjust(hspace=0.0)
        plt.subplots_adjust(right=0.8)
        #plt.tight_layout()
        if lumi_only == True:
            plt.savefig('lumi_time_fits_%s'%(tunnel), dpi = 300, bbox_inches="tight")
            #plt.savefig('test', dpi = 300, bbox_inches="tight")
        plt.show()
        
    def plot_all_luminosity(self, study_period = "Cont_inj", bins = 15):
        if study_period.lower() == "decay":
            bins = 6
        plt.figure(figsize = (14,7.5))
        plt.rc('legend', fontsize=12)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        lumi_data_avg = self.compute_means_and_errs("Lumi", study_period,bins = bins)
        fits, lumi_rates, lumi_rates_err, lumi_rates_scale, fits_corrected = self.measure_and_fit_lumi_bgs(study_period, bins)
        tpcs = ['palila', 'tako', 'elepaio', 'iiwi', 'nene', 'humu']
        pos = ['z = -5.6m', 'z = -8.0m', 'z = -14m', 'z=+6.6m', 'z=+14m', 'z = +16m']
        i=1
        x = np.linspace(0,2,10000) #for plotting fit
        for tpc in tpcs:
            plt.subplot(2,3,i)
            plt.errorbar(lumi_data_avg['ECL_lumi']/10000,lumi_rates[tpc],lumi_rates_err[tpc],lumi_data_avg['ECL_lumi_err']/10000,'o', label = 'Data', alpha = 0.6)
            #plt.plot(x, fits['%s_int'%(tpc)]+ fits['%s_slope'%(tpc)]*x, color = 'tab:blue', label = r'Offset=%s$\pm$%s'%(float('%.2g' % fits[tpc+'_int']), float('%.2g' % fits[tpc+'_int_err'])))
            plt.plot(x, fits['%s_int'%(tpc)]+ fits['%s_slope'%(tpc)]*x, color = 'k', label = r'Fit')
            #plt.fill_between([0],[0],[0], lw = 0, label = r'Slope=%s$\pm$%s'%(float('%.2g' % fits[tpc+'_slope']), float('%.2g' % fits[tpc+'_slope_err'])), color = 'white')
            #if i > 3:
            #plt.errorbar(lumi_data_avg['ECL_lumi']/10000,lumi_rates_scale[tpc],lumi_rates_err[tpc],lumi_data_avg['ECL_lumi_err']/10000,'o', label = 'Corrected', alpha = 0.6)
                #plt.plot(x, fits_corrected['%s_int'%(tpc)]+ fits_corrected['%s_slope'%(tpc)]*x, color = 'tab:orange', label = r'Corrected offset=%s$\pm$%s'%(float('%.2g' % fits_corrected[tpc+'_int']), float('%.2g' % fits_corrected[tpc+'_int_err'])))
                #plt.plot(x, fits['%s_int'%(tpc)]+ fits['%s_slope'%(tpc)]*x, color = 'tab:blue', label = r'Uncorrected fit')
            #plt.plot(x, fits_corrected['%s_int'%(tpc)]+ fits_corrected['%s_slope'%(tpc)]*x, color = 'tab:orange', label = r'Corrected fit')
                #plt.fill_between([0],[0],[0], lw = 0, label = r'Corrected slope=%s$\pm$%s'%(float('%.2g' % fits_corrected[tpc+'_slope']), float('%.2g' % fits_corrected[tpc+'_slope_err'])), color = 'white')
            plt.xlim(0,2)
            plt.xlabel(r'Luminosity [$10^{34}$cm$^{-2}$s$^{-1}$]')
            if i < 4:
                #plt.ylim(-1.5,3.5)
                #plt.yticks([-1,0,1,2,3])
                plt.ylabel(r'R$_L$ [Hz]')
            else:
                #plt.ylim(-0.15,0.35)
                plt.ylabel(r'R$_L$ [Hz]',labelpad = 0.1)
            plt.legend(loc = 'best')
            plt.title(pos[i-1])
            plt.grid()
            i+=1
            plt.subplots_adjust(hspace=0.31)
            plt.subplots_adjust(wspace=0.29)
            plt.subplots_adjust(top=0.935)
            plt.subplots_adjust(bottom=0.090)
            plt.subplots_adjust(left=0.065)
            plt.subplots_adjust(right=0.935)
        plt.tight_layout()
        plt.savefig('lumi_fits_newest.png')
        plt.show()

    def plot_bg_summary(self, study_period, E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0, bins = 6, MC = False, I_HER = 1000, I_LER = 1200, sy_LER=37, sy_HER=36, nb_LER=1576, nb_HER=1576, L=25):
        #tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        plt.rc('legend', fontsize=13)
        plt.rc('xtick', labelsize=20)
        plt.rc('ytick', labelsize=20)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        if MC == True:
            df = self.get_MC_rates(E_cut = E_cut, I_HER = I_HER, I_LER = I_LER, sy_LER=sy_LER, sy_HER=sy_HER, nb_LER=nb_LER, nb_HER=nb_HER, lumi=L)
            #df = self.MC_rates
            df = df[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
            df['total']=df.sum(axis = 1)
            df = df.apply(lambda x: x/df['total'])
            df = df.apply(lambda x: x*100)
            df = df.drop(columns = ['total'])
        else:
            tpcs = ['elepaio', 'tako', 'palila', 'iiwi', 'nene', 'humu']
            LER_fit_params = self.get_fit_parameters("LER", study_period, bins)
            HER_fit_params = self.get_fit_parameters("HER", study_period, bins)
            fit_dict = {}
            df = pd.DataFrame() #order is LER_bg, LER_T, HER_bg, HER_T, Lumi
            lumi_fits = self.measure_and_fit_lumi_bgs(study_period ,bins = 15)[4]
            for tpc in tpcs:
                #I_HER = 1000
                #I_LER = 1200
                #sy_LER = 37
                #sy_HER = 36
                #nb_LER = 1576
                #nb_HER = 1576
                #L = 25 #luminoisty in 1e34 cm-2s-1
                #lumi_fits = self.measure_and_fit_lumi_bgs(tpc, study_period ,bins=20)
                fit_dict['%s_lumi_int'%(tpc)]=lumi_fits['%s_int'%(tpc)]
                fit_dict['%s_lumi_slope'%(tpc)]=lumi_fits['%s_slope'%(tpc)]
                fit_dict['%s_LER_B0'%(tpc)] = LER_fit_params[tpc+'_B0']
                fit_dict['%s_LER_B1'%(tpc)] = LER_fit_params[tpc+'_B1']
                fit_dict['%s_LER_T'%(tpc)] = LER_fit_params[tpc+'_T']
                fit_dict['%s_HER_B0'%(tpc)] = HER_fit_params[tpc+'_B0']
                fit_dict['%s_HER_B1'%(tpc)] = HER_fit_params[tpc+'_B1']
                fit_dict['%s_HER_T'%(tpc)] = HER_fit_params[tpc+'_T']
                LER_bg_base = fit_dict['%s_LER_B0'%(tpc)]*I_LER
                LER_bg_dynamic = fit_dict['%s_LER_B1'%(tpc)]*I_LER**2
                LER_T = fit_dict['%s_LER_T'%(tpc)]*I_LER**2/(sy_LER*nb_LER)
                #print("%s: %s"%(tpc, (LER_bg + LER_T)))
                HER_bg_base = fit_dict['%s_HER_B0'%(tpc)]*I_HER
                HER_bg_dynamic = fit_dict['%s_HER_B1'%(tpc)]*I_HER**2
                HER_T = fit_dict['%s_HER_T'%(tpc)]*I_HER**2/(sy_HER*nb_HER)
                Lumi = fit_dict['%s_lumi_int'%(tpc)]+fit_dict['%s_lumi_slope'%(tpc)]*L
                if Lumi < 0:
                    Lumi = 0
                total = LER_bg_base + LER_bg_dynamic + LER_T + HER_bg_base + HER_bg_dynamic + HER_T + Lumi
                LER_bg_base_frac = LER_bg_base/total * 100
                LER_bg_dynamic_frac = LER_bg_dynamic/total * 100
                HER_bg_base_frac = HER_bg_base/total * 100
                HER_bg_dynamic_frac = HER_bg_dynamic/total * 100
                LER_T_frac = LER_T/total * 100
                HER_T_frac = HER_T/total * 100
                Lumi_frac = Lumi/total * 100
                df[tpc] = [LER_bg_base_frac, LER_bg_dynamic_frac, LER_T_frac, HER_bg_base_frac, HER_bg_dynamic_frac, HER_T_frac, Lumi_frac]
            df = df.T
            df.columns = ['LER Beam Gas Base', 'LER Beam Gas Dyn.', 'LER Touschek', 'HER Beam Gas Base', 'HER Beam Gas Dyn.', 'HER Touschek', 'Luminosity']
            plt.close()
        rc('text', usetex=False)
        plt.rc('legend', fontsize=12)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        colors = ['cyan', 'dodgerblue', 'magenta', 'yellow', 'gold', 'purple', 'limegreen']
        df.index = ['-14m', '-8.0m', '-5.6m', '+6.6m', '+14m', '+16m']
        df.plot(kind='bar', stacked=True, color = colors, legend = False)
        #if study_period == "Cont_inj":
        #    plt.title("Background breakdown continuous injection fits")
        #else:
        #    plt.title("Background breakdown decay fits")
        #plt.xlabel("TPC")
        plt.ylabel('Background Fraction [%]')
        plt.ylim(0,120)
        plt.xticks(rotation=45, ha="right")
        plt.xlabel('TPC z position')
        plt.legend(framealpha = 1, ncol = 4)
        #plt.tight_layout()
        #plt.savefig("bg_breakdown.png")
        plt.show()
        return df

    def compute_data_MC_ratios(self, study_period, E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}, E_cut_err = 0, bins = 6, I_HER = 1000, I_LER = 1200, sy_LER=37, sy_HER=36, nb_LER=1576, nb_HER=1576, L=25):

        MC = self.get_MC_rates(E_cut = E_cut, I_HER = I_HER, I_LER = I_LER, sy_LER=sy_LER, sy_HER=sy_HER, nb_LER=nb_LER, nb_HER=nb_HER, lumi=L)
        #MC = self.MC_rates
        tpcs = ['elepaio', 'tako', 'palila', 'iiwi', 'nene', 'humu']
        LER_fit_params = self.get_fit_parameters("LER", study_period, bins, E_cut=E_cut, E_cut_err = E_cut_err)
        HER_fit_params = self.get_fit_parameters("HER", study_period, bins, E_cut=E_cut, E_cut_err = E_cut_err)
        fit_dict = {}
        df = pd.DataFrame() #order is LER_bg, LER_T, HER_bg, HER_T, Lumi
        lumi_fits = self.measure_and_fit_lumi_bgs(study_period ,bins = 15,E_cut=E_cut, E_cut_err = E_cut_err)[4]
        for tpc in tpcs:
            #I_HER = 1000
            #I_LER = 1200
            #sy_LER = 37
            #sy_HER = 36
            #nb_LER = 1576
            #nb_HER = 1576
            #L = 25 #luminoisty in 1e34 cm-2s-1
            #lumi_fits = self.measure_and_fit_lumi_bgs(tpc, study_period ,bins=20)
            fit_dict['%s_lumi_int'%(tpc)]=lumi_fits['%s_int'%(tpc)]
            fit_dict['%s_lumi_int_err'%(tpc)]=lumi_fits['%s_int_err'%(tpc)]
            fit_dict['%s_lumi_slope'%(tpc)]=lumi_fits['%s_slope'%(tpc)]
            fit_dict['%s_lumi_slope_err'%(tpc)]=lumi_fits['%s_slope_err'%(tpc)]
            fit_dict['%s_LER_B0'%(tpc)] = LER_fit_params[tpc+'_B0']
            fit_dict['%s_LER_B1'%(tpc)] = LER_fit_params[tpc+'_B1']
            fit_dict['%s_LER_T'%(tpc)] = LER_fit_params[tpc+'_T']
            fit_dict['%s_HER_B0'%(tpc)] = HER_fit_params[tpc+'_B0']
            fit_dict['%s_HER_B1'%(tpc)] = HER_fit_params[tpc+'_B1']
            fit_dict['%s_HER_T'%(tpc)] = HER_fit_params[tpc+'_T']
            fit_dict['%s_LER_B0_err'%(tpc)] = LER_fit_params[tpc+'_B0_err']
            fit_dict['%s_LER_B1_err'%(tpc)] = LER_fit_params[tpc+'_B1_err']
            fit_dict['%s_LER_T_err'%(tpc)] = LER_fit_params[tpc+'_T_err']
            fit_dict['%s_HER_B0_err'%(tpc)] = HER_fit_params[tpc+'_B0_err']
            fit_dict['%s_HER_B1_err'%(tpc)] = HER_fit_params[tpc+'_B1_err']
            fit_dict['%s_HER_T_err'%(tpc)] = HER_fit_params[tpc+'_T_err']
            
            LER_bg_base = fit_dict['%s_LER_B0'%(tpc)]*I_LER
            LER_bg_dynamic = fit_dict['%s_LER_B1'%(tpc)]*I_LER**2
            try:
                LER_bg_base_err = LER_bg_base*fit_dict['%s_LER_B0_err'%(tpc)]/fit_dict['%s_LER_B0'%(tpc)]
            except ZeroDivisionError:
                LER_bg_base_err = np.nan
            try:
                LER_bg_dynamic_err = LER_bg_dynamic*fit_dict['%s_LER_B1_err'%(tpc)]/fit_dict['%s_LER_B1'%(tpc)]
            except ZeroDivisionError:
                LER_bg_dynamic_err = np.nan
            HER_bg_base = fit_dict['%s_HER_B0'%(tpc)]*I_HER
            HER_bg_dynamic = fit_dict['%s_HER_B1'%(tpc)]*I_HER**2
            try:
                HER_bg_base_err = HER_bg_base*fit_dict['%s_HER_B0_err'%(tpc)]/fit_dict['%s_HER_B0'%(tpc)]
            except ZeroDivisionError:
                HER_bg_base_err = np.nan
            try:
                HER_bg_dynamic_err = HER_bg_dynamic*fit_dict['%s_HER_B1_err'%(tpc)]/fit_dict['%s_HER_B1'%(tpc)]
            except ZeroDivisionError:
                HER_bg_dynamic_err = np.nan
            
            LER_T = fit_dict['%s_LER_T'%(tpc)]*I_LER**2/(sy_LER*nb_LER)
            try:
                LER_T_err = LER_T*fit_dict['%s_LER_T_err'%(tpc)]/fit_dict['%s_LER_T'%(tpc)]
            except ZeroDivisionError:
                LER_T_err = np.nan
            HER_T = fit_dict['%s_HER_T'%(tpc)]*I_HER**2/(sy_HER*nb_HER)
            try:
                HER_T_err = HER_T*fit_dict['%s_HER_T_err'%(tpc)]/fit_dict['%s_HER_T'%(tpc)]
            except ZeroDivisionError:
                HER_T_err = np.nan
            Lumi = fit_dict['%s_lumi_int'%(tpc)]+fit_dict['%s_lumi_slope'%(tpc)]*L
            if Lumi < 0:
                Lumi = 0
                Lumi_err = 0
            else:
                Lumi_err = np.sqrt((L*fit_dict['%s_lumi_slope_err'%(tpc)])**2 + (fit_dict['%s_lumi_int_err'%(tpc)])**2)
            df[tpc] = [LER_bg_base, LER_bg_base_err, LER_bg_dynamic, LER_bg_dynamic_err, LER_T, LER_T_err, HER_bg_base, HER_bg_base_err, HER_bg_dynamic, HER_bg_dynamic_err, HER_T, HER_T_err, Lumi, Lumi_err]
        df = df.T
        df.columns = ['LER_bg_base',  'LER_bg_base_err', 'LER_bg_dynamic', 'LER_bg_dynamic_err', 'LER_T', 'LER_T_err', 'HER_bg_base', 'HER_bg_base_err', 'HER_bg_dynamic', 'HER_bg_dynamic_err', 'HER_T', 'HER_T_err', 'Lumi', 'Lumi_err']
        plt.close()
        data = df[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
        data_err = df[['LER_bg_base_err', 'LER_bg_dynamic_err', 'LER_T_err', 'HER_bg_base_err', 'HER_bg_dynamic_err', 'HER_T_err', 'Lumi_err']]
        sim = MC[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
        sim_err = MC[['LER_bg_base_err', 'LER_bg_dynamic_err', 'LER_T_err', 'HER_bg_base_err', 'HER_bg_dynamic_err', 'HER_T_err', 'Lumi_err']]
        data_MC = data/sim
        data_MC_err = pd.DataFrame()
        for col in data_MC.columns:
            err = data_MC[col]*np.sqrt((data_err[col+'_err']/data[col])**2+(sim_err[col+'_err']/sim[col])**2)
            data_MC_err[col+'_err'] = err
        data_MC_ratio = pd.concat([data_MC,data_MC_err],axis = 1)
        data_MC_ratio = data_MC_ratio[['LER_bg_base',  'LER_bg_base_err', 'LER_bg_dynamic', 'LER_bg_dynamic_err', 'LER_T', 'LER_T_err', 'HER_bg_base', 'HER_bg_base_err', 'HER_bg_dynamic', 'HER_bg_dynamic_err', 'HER_T', 'HER_T_err', 'Lumi', 'Lumi_err']]
        return df, MC, data_MC_ratio

#E_cut = {'palila': 8.8, 'iiwi': 8.8, 'tako': 5.0, 'nene': 5.6, 'elepaio': 6.0, 'humu': 6.6}
E_cut = {'palila': 8, 'iiwi': 8.8, 'tako': 6, 'nene': 5.5, 'elepaio': 6, 'humu': 6.6}
a = analysis(E_cut = E_cut, E_cut_err = 0)
