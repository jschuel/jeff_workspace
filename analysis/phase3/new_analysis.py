import pandas as pd
import uproot as ur
import root_pandas as rp
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import array
from sklearn import linear_model
from sklearn.linear_model.base import LinearModel
from sklearn.base import RegressorMixin
from sklearn.utils import check_X_y
import numpy as np
import matplotlib
from matplotlib.lines import Line2D
from matplotlib import rc
rc('text', usetex=False)

class analysis:

    def __init__(self, input_file= "~/data/phase3/spring_2020/05-09-20/combined_ntuples/05-09_whole_study_even_newerest.root"):
        self.LER_inj_avg = self.compute_means_and_errs("LER", "Cont_inj")
        self.HER_inj_avg = self.compute_means_and_errs("HER", "Cont_inj")
        self.LER_decay_avg = self.compute_means_and_errs("LER", "Decay")
        self.HER_decay_avg = self.compute_means_and_errs("HER", "Decay")
    
    def get_raw_study_data(self, input_file= "~/data/phase3/spring_2020/05-09-20/combined_ntuples/05-09_whole_study_even_newerest.root"):
        #study_data = rp.read_root(input_file)
        study_data = ur.open(input_file)[ur.open(input_file).keys()[0]].pandas.df(flatten=False)
        return study_data

    def get_tpc_data(self, input_dir = '~/data/phase3/spring_2020/05-09-20/tpc_root_files/'):
        data = {}
        tpcs = ['iiwi', 'humu', 'nene', 'tako', 'palila', 'elepaio']
        for tpc in tpcs:
            #data[tpc] = rp.read_root(input_dir + "%s_all_recoils_only_newest.root"%(tpc))
            data[tpc] = ur.open(input_dir + "%s_all_recoils_only_even_newester.root"%(tpc))[ur.open(input_dir + "%s_all_recoils_only_even_newester.root"%(tpc)).keys()[0]].pandas.df(flatten=False)
        return data

    def get_MC_data(self):
        tpcs = ['iiwi', 'palila', 'tako', 'elepaio']
        bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all', 'RBB_Lumi', 'twoPhoton_Lumi']
        tree = 'tree_fe4_after_threshold'
        data = {}
        truth = {}
        for tpc in tpcs:
            data[tpc] = {}
            truth[tpc] = {}
            dir = '/home/jeef/data/phase3/spring_2020/05-09-20/geant4_simulation/%s/'%(tpc)
            for bg in bgtype:
                try:
                    data[tpc][bg] = ur.open(dir+bg+'_%s.root'%(tpc))[tree].pandas.df(flatten=False)
                    truth[tpc][bg] = ur.open(dir+bg+'_'+tpc+'_'+'truth.root')['recoils'].pandas.df(flatten=False)
                except FileNotFoundError:
                    data[tpc][bg] = pd.DataFrame()
                    truth[tpc][bg] = pd.DataFrame()
                try:
                    truth[tpc][bg] = truth[tpc][bg].loc[truth[tpc][bg].index.isin(data[tpc][bg]['truth_index'])==True]
                    truth[tpc][bg].index = [i for i in range(0,len(truth[tpc][bg]))]
                    data[tpc][bg]['ionization_energy'] = truth[tpc][bg]['eventIonizationEnergy']
                    data[tpc][bg]['truth_energy'] = truth[tpc][bg]['truthRecoilEnergy']
                    data[tpc][bg]['truth_mother_energy'] = truth[tpc][bg]['truthMotherEnergy']
                    data[tpc][bg][['truth_vertex_X', 'truth_vertex_Y', 'truth_vertex_Z', 'truth_px', 'truth_py', 'truth_pz', 'truth_mother_X', 'truth_mother_Y', 'truth_mother_Z', 'truth_mother_px', 'truth_mother_py', 'truth_mother_pz']] = truth[tpc][bg][['truthRecoilVtx_x_belle_frame', 'truthRecoilVtx_y_belle_frame', 'truthRecoilVtx_z_belle_frame', 'truthRecoilMom_x_belle_frame', 'truthRecoilMom_y_belle_frame', 'truthRecoilMom_z_belle_frame', 'truthMotherVtx_x_belle_frame', 'truthMotherVtx_y_belle_frame', 'truthMotherVtx_z_belle_frame', 'truthMotherMom_x_belle_frame', 'truthMotherMom_y_belle_frame', 'truthMotherMom_z_belle_frame']]
                except KeyError:
                    truth[tpc][bg] = pd.DataFrame()
        data_red = {}
        dfs = {}
        for tpc in tpcs:
            data_red[tpc] = {}
            dfs[tpc] = pd.DataFrame()
            for key in data[tpc].keys():
                if len(data[tpc][key]) != 0:
                    data_red[tpc][key] = data[tpc][key]
                    data_red[tpc][key]['bgType'] = key
                    dfs[tpc] = dfs[tpc].append(data_red[tpc][key])
            dfs[tpc].index = [i for i in range(0,len(dfs[tpc]))]
            #dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('Coulomb', 'BeamGas', regex=False)
            #dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('Brems', 'BeamGas', regex=False)
            #dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('_all', '', regex=False)
            #dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('twoPhoton_', '', regex=False)
            #dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('RBB_', '', regex=False)
            
        return dfs

    def apply_energy_calibrations_to_MC(self):
        MC = self.get_MC_data()
        tpcs = MC.keys()
        gain = {'iiwi': 1502, 'palila': 1033, 'tako': 807, 'elepaio': 797}
        W = 35.075
        tot_to_q = {}
        for tpc in tpcs:
            tot_to_q[tpc] = pd.DataFrame()
            tot_to_q[tpc]['tot_code'] = [i for i in range(0,14)]
            
        tot_to_q['iiwi']['conversion'] = np.array([1833.00, 2345.17, 3017.33, 6001.54, 8891.71,
                                  11497.43, 14335.32, 18081.33, 22526.06, 27236.90,
                                  32056.16, 36955.09, 41874.75, 46794.40])

        tot_to_q['palila']['conversion'] = np.array([1768.71, 2202.75, 2670.76, 4049.25, 6586.25,
                                    8954.45, 11551.60, 14428.46, 17618.81, 21140.34,
                                    24831.56, 28804.80, 33534.23, 40821.35])

        tot_to_q['tako']['conversion'] = np.array([2761.20, 3077.66, 3509.80, 5475.02, 9230.59, 
                                  11955.00, 16837.46, 20761.78, 24514.73, 28445.96, 
                                  33071.27, 38033.29, 43011.21, 47989.15])

        tot_to_q['elepaio']['conversion'] = np.array([1859.09, 2496.61, 4128.03, 6844.95, 9450.70,
                                     12158.68, 15125.31, 18507.89, 22166.14, 25826.40,
                                     29597.06, 33588.70, 38207.92, 42827.15])
        for tpc in tpcs:
            MC[tpc]['q_from_tot'] = MC[tpc]['tot']
            
            for i in range(0,len(MC[tpc])):
                try:
                    MC[tpc]['q_from_tot'].iloc[i] = pd.Series(MC[tpc]['tot'].iloc[i]).map(tot_to_q[tpc].set_index('tot_code')['conversion']).to_numpy()
                except ValueError:
                    print(tpc, i)
            MC[tpc]['sumtot'] = [MC[tpc]['q_from_tot'][i].sum() for i in range(0,len(MC[tpc]))]
            MC[tpc]['reco_energy'] = MC[tpc]['sumtot']*35.075/gain[tpc]*1e-3

        return MC

    def get_MC_rates(self):
        tpcs = ['palila', 'tako', 'elepaio', 'iiwi']
        bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all', 'RBB_Lumi', 'twoPhoton_Lumi']
        tree = 'tree_fe4_after_threshold'
        MC = self.apply_energy_calibrations_to_MC()
        rates = {}
        df = pd.DataFrame()
        for tpc in tpcs:
            rates[tpc] = {}
            MC[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>8]
            MC[tpc].index = [i for i in range(0,len(MC[tpc]))]
            for bg in bgtype:
                if (bg == 'Brems_HER_dynamic'):
                    t = 400.
                elif (bg == 'Brems_HER_base'):
                    t = 100.
                elif (bg == 'Coulomb_HER_base') or (bg == 'Coulomb_HER_dynamic'):
                    t = 40.
                elif (bg == 'Brems_LER_base') or (bg == 'Brems_LER_dynamic'):
                    t = 40.
                elif (bg == 'Coulomb_LER_base') or (bg == 'Coulomb_LER_dynamic'):
                    t = 4.0
                elif (bg == 'Touschek_HER_all'):
                    t = 0.8
                elif (bg == 'Touschek_LER_all'):
                    t = 0.1
                elif (bg == 'RBB_Lumi'):
                    if tpc == 'tako':
                        t = 6e-5*(10) # *10 is to accoiunt for reducing the luminosity to levels seen in data
                    elif tpc == 'palila':
                        t = 6e-5*(10) #TODO: Call function that gets lumi slope rather than input by hand
                    elif tpc == 'elepaio':
                        t = 6e-5*(10)
                    elif tpc == 'iiwi':
                        t = 6e-5*(10)
                elif (bg == 'twoPhoton_Lumi'):
                    t = 0.005
                try:
                    rates[tpc][bg] = len(MC[tpc].loc[MC[tpc]['bgType']==bg])/(t*100)
                except FileNotFoundError:
                    rates[tpc][bg] = 0
            df = df.append(pd.DataFrame.from_dict(rates[tpc], 'index').T)
        df.index = tpcs
        df['LER_bg_base'] = df['Brems_LER_base'] + df['Coulomb_LER_base']
        df['HER_bg_base'] = df['Brems_HER_base'] + df['Coulomb_HER_base']
        df['LER_bg_dynamic'] = df['Brems_LER_dynamic'] + df['Coulomb_LER_dynamic']
        df['HER_bg_dynamic'] = df['Brems_HER_dynamic'] + df['Coulomb_HER_dynamic']
        df['LER_T'] = df['Touschek_LER_all']
        df['HER_T'] = df['Touschek_HER_all']
        df['Lumi'] = df['RBB_Lumi'] + df['twoPhoton_Lumi']
        df = df[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
        return df

    def select_study(self, study_type, study_period): #LER, HER, Lumi, Cont_inj, Decay
        raw_data = self.get_raw_study_data()
        study_data = raw_data.loc[(raw_data['%s_study_flag'%(study_type)]==1) & (raw_data['%s_flag'%(study_period)] == 1)]
        return study_data

    def get_tpc_data_during_study_period(self, study_type, study_period):
        study_data = self.select_study(study_type, study_period)
        tpc_data = self.get_tpc_data()
        tpcs = tpc_data.keys()
        tpc_study_data = {}
        for tpc in tpcs:
            tpc_study_data[tpc] = tpc_data[tpc].loc[(tpc_data[tpc]['timestamp_start'].astype('int')).isin(study_data['ts'].astype('int'))==True]
        return tpc_study_data

    def partition_data_into_subsets(self, study_type, study_period, bins=10):
        study_data = self.select_study(study_type, study_period)
        study_data.index = [i for i in range(0,len(study_data))] #reindexes from 0 to length of dataset
        partition_indices = [study_data.index.to_list()[0]] + study_data.loc[np.abs(study_data['ts'].diff())>10].index.to_list() + [study_data.index.to_list()[len(study_data)-1]]
        data_subsets = {}
        for i in range(0,len(partition_indices)-1):
            data_subsets['fill_%s'%(i)] = [i for i in range(partition_indices[i], partition_indices[i+1])]
        #print(data_subsets)
        dfs = {}
        for i in range(0,len(data_subsets)):
            dfs['fill_%s'%(i)] = np.array_split(study_data.iloc[data_subsets['fill_%s'%(i)]], bins)
        return dfs

    def compute_means_and_errs(self, study_type, study_period, bins=10):
        partitioned_data = self.partition_data_into_subsets(study_type, study_period, bins)
        means = partitioned_data['fill_0'][0].head(1) #to give structure to new dataframe
        errs = partitioned_data['fill_0'][0].head(1) #to give structure to new dataframe
        for col in means.columns:
            means[col]=0 #will be deleted at the end. Used so we can append rows of means
            errs[col] = 0
        for key in partitioned_data.keys():
            for i in range(0,len(partitioned_data[key])):
                means = means.append(partitioned_data[key][i].mean(), ignore_index=True)
                errs = errs.append(partitioned_data[key][i].std()/np.sqrt(len(partitioned_data[key][i])), ignore_index=True)
        means = means.loc[means.index>0] # to get rid of phantom row of zeros used for formatting
        means.index = [i for i in range(0,len(means))] #reindex
        errs = errs.loc[errs.index>0] # to get rid of phantom row of zeros used for formatting
        errs.index = [i for i in range(0,len(errs))] #reindex
        errs.columns = [str(col) + '_err' for col in errs.columns]
        for col in errs.columns:
            means[col] = errs[col]
        means = means.drop(columns = ['LER_study_flag_err', 'HER_study_flag_err', 'Lumi_study_flag_err', 'Cont_inj_flag_err', 'Decay_flag_err', 'Nb_HER_err', 'Nb_LER_err'])
        return means

    def get_fit_parameters(self, study_type, study_period, bins=10): #Gives parameters B0, B1, and T defined by Rate/I = B0 + B1*I + T*I/(sy*Nb)
        averaged_data = self.compute_means_and_errs(study_type, study_period, bins)
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
            f2 = ROOT.TF2("f2","[0] + [1]*x + [2]*y", 0, X['I_'+study_type].max(), 0, X['heuristic_x'].max())
            f2.SetParLimits(0,1e-13,1e-5)
            f2.SetParLimits(1,1e-8,1e-3)
            f2.SetParLimits(2,0,1)
            #f2.SetParLimits(0,1e-6,5e-2)
            #f2.SetParLimits(1,1e-15,1)
            #f2.SetParLimits(2,1e-6,1)
            gr = ROOT.TGraph2DErrors(len(x1), x1, x2, y_root, x1err, x2err, y_root_err)
            gr.Fit(f2, 'SREM')
            fit[tpc+'_B0'] = f2.GetParameter(0)
            fit[tpc+'_B1'] = f2.GetParameter(1)
            fit[tpc+'_T'] = f2.GetParameter(2)
            fit[tpc+'_B0_err'] = f2.GetParError(0)
            fit[tpc+'_B1_err'] = f2.GetParError(1)
            fit[tpc+'_T_err'] = f2.GetParError(2)
        return fit

    def plot_fit(self, tpc, study_type,  bins=10, ymax = 2, legend = False):
        study_periods = ["Cont_inj", "Decay"]
        fit_params = {}
        data = {}
        data_avg = {}
        fit_bg_base  = {}
        fit_bg_dynamic  = {}
        fit_bg_base_err = {}
        fit_bg_dynamic_err = {}
        fit_t = {}
        fit_t_err = {}
        fit = {}
        fit_err = {}
        fit_bg_base_avg  = {}
        fit_bg_base_avg_err = {}
        fit_bg_dynamic_avg  = {}
        fit_bg_dynamic_avg_err = {}
        fit_t_avg = {}
        fit_t_avg_err = {}
        fit_avg = {}
        fit_avg_err = {}
        for study_period in study_periods:
            fit_params[study_period] = self.get_fit_parameters(study_type, study_period, bins)
            data[study_period] = self.select_study(study_type, study_period)
            data_avg[study_period] = self.compute_means_and_errs(study_type, study_period, bins)

            
            fit_bg_base[study_period] = fit_params[study_period][tpc+'_B0']*data[study_period]['I_%s'%(study_type)]
            fit_bg_dynamic[study_period] = fit_params[study_period][tpc+'_B1']*data[study_period]['I_%s'%(study_type)]**2
            fit_bg_base_err[study_period] = np.sqrt((fit_params[study_period][tpc+'_B0_err']*data[study_period]['I_%s'%(study_type)])**2)
            fit_bg_dynamic_err[study_period] = np.sqrt((fit_params[study_period][tpc+'_B1_err']*data[study_period]['I_%s'%(study_type)]**2)**2)
            fit_t[study_period] = fit_params[study_period][tpc+'_T']*data[study_period]['I_%s'%(study_type)]**2/(data[study_period]['Sy_%s'%(study_type)]*data[study_period]['Nb_%s'%(study_type)])
            fit_t_err[study_period] = fit_params[study_period][tpc+'_T_err']*data[study_period]['I_%s'%(study_type)]**2/(data[study_period]['Sy_%s'%(study_type)]*data[study_period]['Nb_%s'%(study_type)])
            fit[study_period] = fit_params[study_period][tpc+'_B0']*data[study_period]['I_%s'%(study_type)] + fit_params[study_period][tpc+'_B1']*data[study_period]['I_%s'%(study_type)]**2 + fit_params[study_period][tpc+'_T']*data[study_period]['I_%s'%(study_type)]**2/(data[study_period]['Sy_%s'%(study_type)]*data[study_period]['Nb_%s'%(study_type)])
            fit_err[study_period] = np.sqrt(fit_bg_base_err[study_period]**2 + fit_bg_dynamic_err[study_period]**2 + fit_t_err[study_period]**2)

            
            fit_bg_base_avg[study_period] = fit_params[study_period][tpc+'_B0']*data_avg[study_period]['I_%s'%(study_type)]
            fit_bg_dynamic_avg[study_period] = fit_params[study_period][tpc+'_B1']*data_avg[study_period]['I_%s'%(study_type)]**2
            fit_bg_base_avg_err[study_period] = np.sqrt((fit_params[study_period][tpc+'_B0_err']*data_avg[study_period]['I_%s'%(study_type)])**2)
            fit_bg_dynamic_avg_err[study_period] = np.sqrt((fit_params[study_period][tpc+'_B1_err']*data_avg[study_period]['I_%s'%(study_type)]**2)**2)
            fit_t_avg[study_period] = fit_params[study_period][tpc+'_T']*data_avg[study_period]['I_%s'%(study_type)]**2/(data_avg[study_period]['Sy_%s'%(study_type)]*data_avg[study_period]['Nb_%s'%(study_type)])
            fit_t_avg_err[study_period] = fit_params[study_period][tpc+'_T_err']*data_avg[study_period]['I_%s'%(study_type)]**2/(data_avg[study_period]['Sy_%s'%(study_type)]*data_avg[study_period]['Nb_%s'%(study_type)])
            fit_avg[study_period] = fit_params[study_period][tpc+'_B0']*data_avg[study_period]['I_%s'%(study_type)] + fit_params[study_period][tpc+'_B1']*data_avg[study_period]['I_%s'%(study_type)]**2 + fit_params[study_period][tpc+'_T']*data_avg[study_period]['I_%s'%(study_type)]**2/(data_avg[study_period]['Sy_%s'%(study_type)]*data_avg[study_period]['Nb_%s'%(study_type)])
            fit_avg_err[study_period] = np.sqrt(fit_bg_base_avg_err[study_period]**2 + fit_bg_dynamic_avg_err[study_period]**2 + fit_t_avg_err[study_period]**2)

            
            if study_type == "LER":
                plt.plot(data[study_period]['ts'], data[study_period]['I_%s'%(study_type)], 'o', markersize = 1, color = 'red', label = "I_LER [mA]")
            else:
                plt.plot(data[study_period]['ts'], data[study_period]['I_%s'%(study_type)], 'o', markersize = 1, color = 'blue', label = "I_HER [mA]")
        plt.ylim(0,620)
        plt.ylabel('Beam Current [mA]')
        plt.twinx()
        for study_period in study_periods:
            p1 = plt.errorbar(data_avg[study_period]['ts'], data_avg[study_period][tpc+'_neutrons'], data_avg[study_period][tpc+'_neutrons_err'], data_avg[study_period]['ts_err'], 'o', markersize = 2, color = 'black', label = 'data')
            p2 = plt.errorbar(data_avg[study_period]['ts'], fit_avg[study_period], fit_avg_err[study_period], data_avg[study_period]['ts_err'], 's', markersize = 2, color = 'cyan', label = 'total fit', alpha = 0.6)
            p3 = plt.errorbar(data_avg[study_period]['ts'], fit_bg_base_avg[study_period], fit_bg_base_avg_err[study_period], data_avg[study_period]['ts_err'], 's', markersize = 2, color = 'purple', label = 'beam gas base', alpha = 0.6)
            p4 = plt.errorbar(data_avg[study_period]['ts'], fit_bg_dynamic_avg[study_period], fit_bg_dynamic_avg_err[study_period], data_avg[study_period]['ts_err'], 's', markersize = 2, color = 'gray', label = 'beam gas dyn.', alpha = 0.6)
            p5 = plt.errorbar(data_avg[study_period]['ts'], fit_t_avg[study_period], fit_t_avg_err[study_period], data_avg[study_period]['ts_err'], 's', markersize = 2, color = 'green', label = 'Touschek', alpha = 0.6)

            #p1 = plt.errorbar(data_avg[study_period]['ts'], data_avg[study_period][tpc+'_neutrons'], data_avg[study_period][tpc+'_neutrons_err'], data_avg[study_period]['ts_err'], 'o', markersize = 2, color = 'black', label = 'data')
            #p2 = plt.errorbar(data[study_period]['ts'], fit[study_period], fit_err[study_period], [0 for i in range(0,len(data[study_period]))], 's', markersize = 2, color = 'cyan', label = 'total fit', alpha = 0.3)
            #p3 = plt.errorbar(data[study_period]['ts'], fit_bg_base[study_period], fit_bg_base_err[study_period], [0 for i in range(0,len(data[study_period]))], 's', markersize = 2, color = 'purple', label = 'beam gas base', alpha = 0.3)
            #p4 = plt.errorbar(data[study_period]['ts'], fit_bg_dynamic[study_period], fit_bg_dynamic_err[study_period], [0 for i in range(0,len(data[study_period]))], 's', markersize = 2, color = 'gray', label = 'beam gas dyn.', alpha = 0.3)
            #p5 = plt.errorbar(data[study_period]['ts'], fit_t[study_period], fit_t_err[study_period], [0 for i in range(0,len(data[study_period]))], 's', markersize = 2, color = 'green', label = 'Touschek', alpha = 0.3)
            
            #partitions = [data[study_period].index.to_list()[0]] + data[study_period].loc[np.abs(data[study_period]['ts'].diff())>100].index.to_list() + [data[study_period].index.to_list()[len(data[study_period])-1]]
            #indices = []
            #for i in range(0,len(partitions)-1):
            #    indices = [i for i in range(partitions[i], partitions[i+1])]
            #    data_tmp = data[study_period].loc[data[study_period].index.isin(indices) == True]
            #    fit_bg_tmp = fit_bg[study_period].loc[fit_bg[study_period].index.isin(indices) == True]
            #    fit_bg_err_tmp = fit_bg_err[study_period].loc[fit_bg_err[study_period].index.isin(indices) == True]
            #    fit_t_tmp = fit_t[study_period].loc[fit_t[study_period].index.isin(indices) == True]
            #    fit_t_err_tmp = fit_t_err[study_period].loc[fit_t_err[study_period].index.isin(indices) == True]
            #    fit_tmp = fit[study_period].loc[fit[study_period].index.isin(indices) == True]
            #    fit_err_tmp = fit_err[study_period].loc[fit_err[study_period].index.isin(indices) == True]
            #    plt.fill_between(data_tmp['ts'], (fit_bg_tmp-fit_bg_err_tmp), (fit_bg_tmp+fit_bg_err_tmp), 'o', color = 'cyan', alpha = 0.3, label = 'Predicted beam-gas rate')
            #   plt.fill_between(data_tmp['ts'], (fit_t_tmp-fit_t_err_tmp), (fit_t_tmp+fit_t_err_tmp), 'o', color = 'limegreen', alpha = 0.3, label = 'Predicted Touschek rate')
            #    plt.fill_between(data_tmp['ts'], (fit_tmp-fit_err_tmp), (fit_tmp+fit_err_tmp), 'o', color = 'magenta', alpha = 0.3, label = 'Predicted total rate')
        plt.ylim(-0.01,ymax)
        plt.ylabel('Rate[Hz]')
        if legend == True:
            #if study_type == "LER":
            #    custom_lines = [Line2D([0], [0], color='red', lw=4, label = '%s Current [mA]'%(study_type)), Line2D([0], [0], color='magenta', lw=4, label = 'TPC predicted rate'), Line2D([0], [0], color='limegreen', lw=4, label = 'TPC predicted Touschek rate'), Line2D([0], [0], color='cyan', lw=4, label = 'TPC predicted beam-gas rate'), Line2D([0], [0], marker='o', color='w', markersize=6, lw=0, markerfacecolor = 'black', label = 'Data')]
            #else:
            #    custom_lines = [Line2D([0], [0], color='blue', lw=4, label = '%s Current [mA]'%(study_type)), Line2D([0], [0], color='magenta', lw=4, label = 'TPC predicted rate'), Line2D([0], [0], color='limegreen', lw=4, label = 'TPC predicted Touschek rate'), Line2D([0], [0], color='cyan', lw=4, label = 'TPC predicted beam-gas rate'), Line2D([0], [0], marker='o', color='w', markersize=6, lw=0, markerfacecolor = 'black', label = 'Data')]
            #plt.legend(handles = custom_lines, loc = 'upper right')
            plt.legend([p1,p2,p3,p4,p5],['Data', 'Total Fit', 'Beam Gas Base', 'Beam Gas Dyn.', 'Touschek'])
    def make_fwd_plots(self, study_type, bins=10):
        if study_type == "LER":
            ymax = 2
        else:
            ymax = 2
        plt.figure(figsize = (12,8))
        plt.rc('legend', fontsize=10)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        #plt.subplot(3,1,1)
        self.plot_fit("iiwi", study_type, bins, ymax=1, legend = True)
        plt.title("Iiwi (z = +6.6 m)")
        plt.xticks([])
        #plt.subplot(3,1,2)
        #self.plot_fit("nene", study_type, bins, ymax)
        #plt.title("Nene (z = +14 m)")
        #plt.xticks([])
        #plt.subplot(3,1,3)
        #self.plot_fit("humu", study_type, bins, ymax)
        #plt.title("Humu (z = +16 m)")
        #plt.xticks([])
        #plt.tight_layout()
        #plt.subplots_adjust(hspace = 0.2)
        #plt.subplots_adjust(wspace = 0.5)
        plt.show()

    def make_bwd_plots(self, study_type, bins=10):
        if study_type == "HER":
            ymax = 0.3
        else:
            ymax = 0.3
        plt.figure(figsize = (16,8))
        plt.rc('legend', fontsize=10)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        plt.subplot(3,1,1)
        self.plot_fit("palila", study_type, bins, ymax=0.2, legend = True)
        plt.title("Palila (z = -5.6 m)")
        plt.xticks([])
        plt.subplot(3,1,2)
        self.plot_fit("tako", study_type, bins, ymax=0.2)
        plt.title("Tako (z = -8.0 m)")
        plt.xticks([])
        plt.subplot(3,1,3)
        self.plot_fit("elepaio", study_type, bins, ymax=0.6)
        plt.title("elepaio (z = -14 m)")
        plt.xticks([])
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0.2)
        plt.show()

    def measure_and_fit_lumi_bgs(self, tpc, study_period ,bins=10):
        lumi_data_avg = self.compute_means_and_errs("Lumi", study_period, bins=20)
        LER_fit_params = self.get_fit_parameters("LER", study_period, bins)
        HER_fit_params = self.get_fit_parameters("HER", study_period, bins)
        LER_rates = LER_fit_params[tpc+'_B0']*lumi_data_avg['I_LER'] + LER_fit_params[tpc+'_B1']*lumi_data_avg['I_LER']**2 + LER_fit_params[tpc+'_T']*lumi_data_avg['I_LER']**2/(lumi_data_avg['Sy_LER']*lumi_data_avg['Nb_LER'])
        HER_rates = HER_fit_params[tpc+'_B0']*lumi_data_avg['I_HER'] + HER_fit_params[tpc+'_B1']*lumi_data_avg['I_HER']**2 + HER_fit_params[tpc+'_T']*lumi_data_avg['I_HER']**2/(lumi_data_avg['Sy_HER']*lumi_data_avg['Nb_HER'])
        LER_rates_err = np.sqrt((LER_fit_params[tpc+'_B0']+2*LER_fit_params[tpc+'_B1']*lumi_data_avg['I_LER']+2*LER_fit_params[tpc+'_T']*lumi_data_avg['I_LER']/(lumi_data_avg['Sy_LER']**2*lumi_data_avg['Nb_LER'])*lumi_data_avg['I_LER_err'])**2+(LER_fit_params[tpc+'_T']*lumi_data_avg['I_LER']**2/(lumi_data_avg['Sy_LER']**2*lumi_data_avg['Nb_LER'])*lumi_data_avg['Sy_LER_err'])**2)
        HER_rates_err = np.sqrt((HER_fit_params[tpc+'_B0']+2*HER_fit_params[tpc+'_B1']*lumi_data_avg['I_HER']+2*HER_fit_params[tpc+'_T']*lumi_data_avg['I_HER']/(lumi_data_avg['Sy_HER']**2*lumi_data_avg['Nb_HER'])*lumi_data_avg['I_HER_err'])**2+(HER_fit_params[tpc+'_T']*lumi_data_avg['I_HER']**2/(lumi_data_avg['Sy_HER']**2*lumi_data_avg['Nb_HER'])*lumi_data_avg['Sy_HER_err'])**2)
        lumi_rates = lumi_data_avg[tpc+'_neutrons'] - LER_rates - HER_rates
        lumi_rates_err = np.sqrt(lumi_data_avg[tpc+'_neutrons_err']**2 + LER_rates_err**2 + HER_rates_err**2)
        
        gr = ROOT.TGraphErrors(len(lumi_rates), array.array('d', lumi_data_avg['ECL_lumi']/10000), array.array('d', lumi_rates), array.array('d', lumi_data_avg['ECL_lumi_err']/10000), array.array('d', lumi_rates_err))
        f1 = ROOT.TF1("f1", "[0] + [1]*x", 0, 2)
        gr.Fit("f1", "SEMR")
        fits = {}
        fits['%s_int'%(tpc)] = gr.GetFunction("f1").GetParameter(0)
        fits['%s_int_err'%(tpc)] = gr.GetFunction("f1").GetParError(0)
        fits['%s_slope'%(tpc)] = gr.GetFunction("f1").GetParameter(1)
        fits['%s_slope_err'%(tpc)] = gr.GetFunction("f1").GetParError(1)

        x = np.linspace(0,2,10000) #for plotting fit
        plt.errorbar(lumi_data_avg['ECL_lumi']/10000,lumi_rates,lumi_rates_err,lumi_data_avg['ECL_lumi_err']/10000,'o', label = 'Data', alpha = 0.6)
        plt.plot(x, fits['%s_int'%(tpc)]+ fits['%s_slope'%(tpc)]*x, color = 'black', label = r'offset=%s$\pm$%s'%(float('%.2g' % fits[tpc+'_int']), float('%.2g' % fits[tpc+'_int_err'])))
        plt.fill_between([0],[0],[0], lw = 0, label = r'slope=%s$\pm$%s'%(float('%.2g' % fits[tpc+'_slope']), float('%.2g' % fits[tpc+'_slope_err'])), color = 'white')
        plt.xlim(0,2)
        plt.xlabel(r'Luminosity [$10^{34}$cm$^{-2}$s$^{-1}$]')
        #plt.ylabel(r'Rate - $LER_{fit}$ - $HER_{fit}$ [Hz]')
        plt.ylabel(r'R$_L$ [Hz]')
        plt.ylim(-1.5,3.5)
        plt.yticks([-1,0,1,2,3])
        plt.legend(loc = 'best')
        return fits

    def plot_all_luminosity(self, study_period = "Cont_inj", bins = 20):
        plt.figure(figsize = (15,10))
        plt.rc('legend', fontsize=12)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        plt.subplot(3,2,1)
        self.measure_and_fit_lumi_bgs("palila", study_period, bins)
        plt.title("Palila (z = -5.6 m)")
        plt.grid()
        plt.subplot(3,2,2)
        self.measure_and_fit_lumi_bgs("iiwi", study_period, bins)
        plt.title("Iiwi (z = +6.6 m)")
        plt.grid()
        plt.subplot(3,2,3)
        self.measure_and_fit_lumi_bgs("tako", study_period, bins)
        plt.title("Tako (z = -8.0 m)")
        plt.grid()
        #plt.subplot(3,2,4)
        #self.measure_and_fit_lumi_bgs("nene", study_period, bins)
        #plt.title("Nene (z = 14 m)")
        #plt.grid()
        plt.subplot(3,2,5)
        self.measure_and_fit_lumi_bgs("elepaio", study_period, bins)
        plt.title("Elepaio (z = -14 m)")
        plt.grid()
        #plt.subplot(3,2,6)
        #self.measure_and_fit_lumi_bgs("humu", study_period, bins)
        #plt.title("Humu (z = +16 m)")
        plt.subplots_adjust(hspace=0.5)
        plt.subplots_adjust(wspace=0.5)
        #plt.grid()
        plt.savefig('lumi_fits.png', bbox_inches='tight')
        plt.show()

    def plot_bwd_luminosity(self, bins=10):
        plt.figure(figsize = (16,8))
        plt.rc('legend', fontsize=12)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        plt.subplot(2,3,1)
        self.measure_and_fit_lumi_bgs("palila", "Cont_inj", bins)
        plt.title("Palila (z = -5.6 m) Cont. Inj. Fit")
        plt.subplot(2,3,4)
        self.measure_and_fit_lumi_bgs("palila", "Decay", bins)
        plt.title("Palila (z = -5.6 m) Decay Fit")
        plt.subplot(2,3,2)
        self.measure_and_fit_lumi_bgs("tako", "Cont_inj", bins)
        plt.title("Tako (z = -8.0 m) Cont. Inj. Fit")
        plt.subplot(2,3,5)
        self.measure_and_fit_lumi_bgs("tako", "Decay", bins)
        plt.title("Tako (z = -8.0 m) Decay Fit")
        plt.subplot(2,3,3)
        self.measure_and_fit_lumi_bgs("elepaio", "Cont_inj", bins)
        plt.title("Elepaio (z = -14 m) Cont. Inj. Fit")
        plt.subplot(2,3,6)
        self.measure_and_fit_lumi_bgs("elepaio", "Decay", bins)
        plt.title("Elepaio (z = -14 m) Decay Fit")
        plt.tight_layout()
        plt.show()

    def plot_fwd_luminosity(self, bins=10):
        plt.figure(figsize = (16,8))
        plt.rc('legend', fontsize=12)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        plt.subplot(2,3,1)
        self.measure_and_fit_lumi_bgs("iiwi", "Cont_inj", bins)
        plt.title("Iiwi (z = +6.6 m) Cont. Inj. Fit")
        plt.subplot(2,3,4)
        self.measure_and_fit_lumi_bgs("iiwi", "Decay", bins)
        plt.title("Iiwi (z = +6.6 m) Decay Fit")
        plt.subplot(2,3,2)
        self.measure_and_fit_lumi_bgs("nene", "Cont_inj", bins)
        plt.title("Nene (z = +14 m) Cont. Inj. Fit")
        plt.subplot(2,3,5)
        self.measure_and_fit_lumi_bgs("nene", "Decay", bins)
        plt.title("Nene (z = +14 m) Decay Fit")
        plt.subplot(2,3,3)
        self.measure_and_fit_lumi_bgs("humu", "Cont_inj", bins)
        plt.title("Humu (z = +16 m) Cont. Inj. Fit")
        plt.subplot(2,3,6)
        self.measure_and_fit_lumi_bgs("humu", "Decay", bins)
        plt.title("Humu (z = +16 m) Decay Fit")
        plt.tight_layout()
        plt.show()

    def plot_bg_summary(self, study_period, bins=10, MC = False):
        #tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        if MC == True:
            df = self.get_MC_rates()
            df['total']=df.sum(axis = 1)
            df = df.apply(lambda x: x/df['total'])
            df = df.apply(lambda x: x*100)
            df = df.drop(columns = ['total'])
        else:
            tpcs = ['palila', 'tako', 'elepaio', 'iiwi']
            LER_fit_params = self.get_fit_parameters("LER", study_period, bins)
            HER_fit_params = self.get_fit_parameters("HER", study_period, bins)
            fit_dict = {}
            df = pd.DataFrame() #order is LER_bg, LER_T, HER_bg, HER_T, Lumi 
            for tpc in tpcs:
                I_HER = 1000
                I_LER = 1200
                sy_LER = 37
                sy_HER = 36
                nb_LER = 1576
                nb_HER = 1576
                L = 2.5 #luminoisty in 1e34 cm-2s-1
                lumi_fits = self.measure_and_fit_lumi_bgs(tpc, study_period ,bins=10)
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
        df.plot(kind='bar', stacked=True, color = colors, legend = False)
        #if study_period == "Cont_inj":
        #    plt.title("Background breakdown continuous injection fits")
        #else:
        #    plt.title("Background breakdown decay fits")
        #plt.xlabel("TPC")
        plt.ylabel('Background Fraction [%]')
        plt.ylim(0,120)
        plt.xticks(rotation=45, ha="right")
        plt.legend(framealpha = 1, ncol = 4)
        #plt.tight_layout()
        #plt.savefig("bg_breakdown.png")
        plt.show()
        return df

    def compute_data_MC_ratios(self, study_period, bins=10):

        MC = self.get_MC_rates()

        tpcs = ['palila', 'tako', 'elepaio', 'iiwi']
        LER_fit_params = self.get_fit_parameters("LER", study_period, bins)
        HER_fit_params = self.get_fit_parameters("HER", study_period, bins)
        fit_dict = {}
        df = pd.DataFrame() #order is LER_bg, LER_T, HER_bg, HER_T, Lumi 
        for tpc in tpcs:
            I_HER = 1000
            I_LER = 1200
            sy_LER = 37
            sy_HER = 36
            nb_LER = 1576
            nb_HER = 1576
            L = 2.5 #luminoisty in 1e34 cm-2s-1
            lumi_fits = self.measure_and_fit_lumi_bgs(tpc, study_period ,bins=10)
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
            LER_bg_base_err = LER_bg_base*fit_dict['%s_LER_B0_err'%(tpc)]/fit_dict['%s_LER_B0'%(tpc)]
            LER_bg_dynamic_err = LER_bg_dynamic*fit_dict['%s_LER_B1_err'%(tpc)]/fit_dict['%s_LER_B1'%(tpc)]
            HER_bg_base = fit_dict['%s_HER_B0'%(tpc)]*I_HER
            HER_bg_dynamic = fit_dict['%s_HER_B1'%(tpc)]*I_HER**2
            HER_bg_base_err = HER_bg_base*fit_dict['%s_HER_B0_err'%(tpc)]/fit_dict['%s_HER_B0'%(tpc)]
            HER_bg_dynamic_err = HER_bg_dynamic*fit_dict['%s_HER_B1_err'%(tpc)]/fit_dict['%s_HER_B1'%(tpc)]
            
            LER_T = fit_dict['%s_LER_T'%(tpc)]*I_LER**2/(sy_LER*nb_LER)
            LER_T_err = LER_T*fit_dict['%s_LER_T_err'%(tpc)]/fit_dict['%s_LER_T'%(tpc)]
            HER_T = fit_dict['%s_HER_T'%(tpc)]*I_HER**2/(sy_HER*nb_HER)
            HER_T_err = HER_T*fit_dict['%s_HER_T_err'%(tpc)]/fit_dict['%s_HER_T'%(tpc)]
            
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
        return df, MC
        
    '''
    def get_fit_parameters(self, study_type, study_period, bins=n):
        data = self.heuristic_averages(study_type, study_period, bins)
        tpcs = ['iiwi', 'humu', 'nene', 'tako', 'elepaio', 'palila']
        #x = {}
        #x_err = {}
        #y = {}
        #y_err = {}
        #gr = {}
        #fits = {}
        #for tpc in tpcs:
        #    index = data.loc[data['%s_heuristic_y'%(tpc)]!=0].index.to_numpy()
        #    x[tpc] = array.array('d', data['heuristic_x'][index].to_numpy())
        #    x_err[tpc] = array.array('d', data['heuristic_x_err'][index].to_numpy())
        #    y[tpc] = array.array('d',data['%s_heuristic_y'%(tpc)][index])
        #    y_err[tpc] = array.array('d', data['%s_heuristic_y_err'%(tpc)][index])
        #    gr[tpc] = ROOT.TGraph(len(x[tpc]), x[tpc], y[tpc])#, x_err[tpc], y_err[tpc])
        #    f1 = ROOT.TF1("f1", "[0] + [1]*x", 0,  data['I_'+study_type].max())
            #f1.SetParLimits(1,0,0.5)
            #f1.SetParLimits(0,0,1e5)
        #    gr[tpc].Fit("f1", "SEM")
        #    fits['%s_B'%(tpc)] = gr[tpc].GetFunction("f1").GetParameter(0)
        #    fits['%s_B_err'%(tpc)] = gr[tpc].GetFunction("f1").GetParError(0)
        #    fits['%s_T'%(tpc)] = gr[tpc].GetFunction("f1").GetParameter(1)
        #    fits['%s_T_err'%(tpc)] = gr[tpc].GetFunction("f1").GetParError(1)
        #return fits

    def plot_heuristic_fit(self, tpc, study_type, study_period, bins=n):
        data = self.heuristic_averages(study_type, study_period, bins)
        #fit = self.get_fit_parameters(study_type, study_period, bins)
        Nbs = [393, 783, 1565]
        index_393 = data.loc[data['Nb_'+study_type]==393].index.to_numpy()
        index_783 = data.loc[data['Nb_'+study_type]==783].index.to_numpy()
        index_1565 = data.loc[data['Nb_'+study_type]==1565].index.to_numpy()
        plt.errorbar(data['heuristic_x'][index_393], data[tpc+'_heuristic_y'][index_393], data[tpc+'_heuristic_y_err'][index_393], data['heuristic_x_err'][index_393], 'o', markersize = 3, color = 'magenta', label = 'Nb 393')
        plt.errorbar(data['heuristic_x'][index_783], data[tpc+'_heuristic_y'][index_783], data[tpc+'_heuristic_y_err'][index_783], data['heuristic_x_err'][index_783], 'o', markersize = 3, color = 'tab:blue', label = 'Nb 783')
        plt.errorbar(data['heuristic_x'][index_1565], data[tpc+'_heuristic_y'][index_1565], data[tpc+'_heuristic_y_err'][index_1565], data['heuristic_x_err'][index_1565], 'o', markersize = 3, color = 'tab:green', label = 'Nb 1565')
        plt.plot(data['heuristic_x'],data['heuristic_x']*fit[tpc+'_T']+fit[tpc+'_B'], color = 'black', label = 'B=%s+/-%s, T=%s+/-%s'%(float('%.2g' % fit[tpc+'_B']), float('%.2g' % fit[tpc+'_B_err']), float('%.2g' % fit[tpc+'_T']), float('%.2g' % fit[tpc+'_T_err'])))
        plt.legend()

    def make_fwd_plots(self, study_type, bins=10):
        plt.subplot(2,3,1)
        self.plot_heuristic_fit("iiwi", study_type, "Decay", bins)
        plt.title("Iiwi (z = +6.6m) Decay")
        #plt.xlim(0,800000)
        plt.subplot(2,3,4)
        self.plot_heuristic_fit("iiwi", study_type, "Cont_inj", bins)
        plt.title("Iiwi (z = +6.6m) Cont. Inj.")
        #plt.xlim(0,800000)
        plt.subplot(2,3,2)
        self.plot_heuristic_fit("nene", study_type, "Decay", bins)
        plt.title("Nene (z = +14m) Decay")
        #plt.xlim(0,800000)
        plt.subplot(2,3,5)
        self.plot_heuristic_fit("nene", study_type, "Cont_inj", bins)
        plt.title("Nene (z = +14m) Cont. Inj.")
        #plt.xlim(0,800000)
        plt.subplot(2,3,3)
        self.plot_heuristic_fit("humu", study_type, "Decay", bins)
        plt.title("Humu (z = +16m) Decay")
        #plt.xlim(0,800000)
        plt.subplot(2,3,6)
        self.plot_heuristic_fit("humu", study_type, "Cont_inj", bins)
        plt.title("Humu (z = +16m) Cont. Inj.")
        #plt.xlim(0,800000)
        plt.show()

    def make_bwd_plots(self, study_type, bins=10):
        plt.subplot(2,3,1)
        self.plot_heuristic_fit("palila", study_type, "Decay", bins)
        plt.title("palila (z = -5.6m) Decay")
        #plt.xlim(0,800000)
        plt.subplot(2,3,4)
        self.plot_heuristic_fit("palila", study_type, "Cont_inj", bins)
        plt.title("Palila (z = -5.6m) Cont. Inj.")
        #plt.xlim(0,800000)
        plt.subplot(2,3,2)
        self.plot_heuristic_fit("tako", study_type, "Decay", bins)
        plt.title("Tako (z = -8.0m) Decay")
        #plt.xlim(0,800000)
        plt.subplot(2,3,5)
        self.plot_heuristic_fit("tako", study_type, "Cont_inj", bins)
        plt.title("Tako (z = -8.0m) Cont. Inj.")
        #plt.xlim(0,800000)
        plt.subplot(2,3,3)
        self.plot_heuristic_fit("elepaio", study_type, "Decay", bins)
        plt.title("Elepaio (z = +14m) Decay")
        #plt.xlim(0,800000)
        plt.subplot(2,3,6)
        self.plot_heuristic_fit("elepaio", study_type, "Cont_inj", bins)
        plt.title("Elepaio (z = -14m) Cont. Inj.")
        #plt.xlim(0,800000)
        plt.show()
        
    def plot_fits(self, study_type, study_period, bins=10):
        raw_data = self.compute_means_and_errs(study_type, study_period)
        fits = self.get_fit_parameters(study_type, study_period, bins)
        #plt.plot(raw_data['ts'], raw_data['I_'+study_type], 'o', markersize = 1)
        #plt.errorbar(raw_data['ts'], raw_data['elepaio_neutrons'], raw_data['elepaio_neutrons_err'], 'o', markersize = 1)
        plt.plot(raw_data['I_'+study_type], raw_data['nene_neutrons'], 'o')
    '''
    
a = analysis()

