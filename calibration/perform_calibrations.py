### Program has a simulation class and a calibration class. Write longer description later

import numpy as np
import pandas as pd
import uproot as ur
import root_pandas as rp
import ROOT
from ROOT import TVector3
import array
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc
from matplotlib.patches import Patch
rc('text', usetex=False)

pd.set_option('mode.chained_assignment', None) #suppress error

class simulation:

    def __init__(self):
        #self.make_MC_summary_plots(recoil_species = 'CO', zmin = 0, zmax = 10)
        #self.make_MC_summary_plots(recoil_species = 'He', zmin = 0, zmax = 10)
        #self.make_MC_summary_plots(recoil_species = 'CO', zmin = 2, zmax = 8)
        #self.make_MC_summary_plots(recoil_species = 'He', zmin = 2, zmax = 8)
        pass
    def load_MC(self, base_path = '~/data/phase3/simulation/resolution_paper/tpc_sims/', recoil_species = 'all'):
        after_threshold_C = {}
        after_threshold_O = {}
        after_threshold_He = {}
        after_threshold = {}
        after_drift_C = {}
        after_drift_O = {}
        after_drift_He = {}
        after_drift = {}

        tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        for tpc in tpcs:
            after_threshold_C[tpc] = ur.open(base_path + "%s_after_threshold_C_all.root"%(tpc))[ur.open(base_path + "%s_after_threshold_C_all.root"%(tpc)).keys()[0]].pandas.df(flatten=False)
            after_threshold_O[tpc] = ur.open(base_path + "%s_after_threshold_O_all.root"%(tpc))[ur.open(base_path + "%s_after_threshold_O_all.root"%(tpc)).keys()[0]].pandas.df(flatten=False)
            after_threshold_He[tpc] = ur.open(base_path + "%s_after_threshold_He_all.root"%(tpc))[ur.open(base_path + "%s_after_threshold_He_all.root"%(tpc)).keys()[0]].pandas.df(flatten=False)
            after_drift_C[tpc] = ur.open(base_path + "%s_after_drift_C_all.root"%(tpc))[ur.open(base_path + "%s_after_drift_C_all.root"%(tpc)).keys()[0]].pandas.df(flatten=False)
            after_drift_O[tpc] = ur.open(base_path + "%s_after_drift_O_all.root"%(tpc))[ur.open(base_path + "%s_after_drift_O_all.root"%(tpc)).keys()[0]].pandas.df(flatten=False)
            after_drift_He[tpc] = ur.open(base_path + "%s_after_drift_He_all.root"%(tpc))[ur.open(base_path + "%s_after_drift_He_all.root"%(tpc)).keys()[0]].pandas.df(flatten=False)
            after_threshold_C[tpc]['recoil_species'] = 'C'
            after_threshold_O[tpc]['recoil_species'] = 'O'
            after_threshold_He[tpc]['recoil_species'] = 'He'
            after_drift_C[tpc]['recoil_species'] = 'C'
            after_drift_O[tpc]['recoil_species'] = 'O'
            after_drift_He[tpc]['recoil_species'] = 'He'
            after_threshold[tpc] = after_threshold_He[tpc].append(after_threshold_C[tpc].append(after_threshold_O[tpc]))
            after_threshold[tpc] = after_threshold[tpc].sort_values(by = 'truth_energy')
            after_drift[tpc] = after_drift_He[tpc].append(after_drift_C[tpc].append(after_drift_O[tpc]))
            after_drift[tpc] = after_drift[tpc].sort_values(by = 'truth_energy')
            if recoil_species == 'He':
                after_threshold[tpc] = after_threshold[tpc].loc[after_threshold[tpc]['recoil_species'] == 'He']
                after_drift[tpc] = after_drift[tpc].loc[after_drift[tpc]['recoil_species'] == 'He']
            elif recoil_species == 'C':
                after_threshold[tpc] = after_threshold[tpc].loc[after_threshold[tpc]['recoil_species'] == 'C']
                after_drift[tpc] = after_drift[tpc].loc[after_drift[tpc]['recoil_species'] == 'C']
            elif recoil_species == 'O':
                after_threshold[tpc] = after_threshold[tpc].loc[after_threshold[tpc]['recoil_species'] == 'O']
                after_drift[tpc] = after_drift[tpc].loc[after_drift[tpc]['recoil_species'] == 'O']
            elif recoil_species == 'CO':
                after_threshold[tpc] = after_threshold[tpc].loc[(after_threshold[tpc]['recoil_species'] == 'O')
                                                                | (after_threshold[tpc]['recoil_species'] == 'C')]
                after_drift[tpc] = after_drift[tpc].loc[(after_drift[tpc]['recoil_species'] == 'O')
                                                        | (after_drift[tpc]['recoil_species'] == 'C')]
            after_threshold[tpc].index = [i for i in range(0,len(after_threshold[tpc]))] #reindex
            after_drift[tpc].index = [i for i in range(0,len(after_drift[tpc]))] #reindex 
        return after_threshold, after_drift

    def apply_energy_calibrations(self, base_path =  '~/data/phase3/simulation/resolution_paper/tpc_sims/', recoil_species = 'all', zmin = 0, zmax = 10):
        MC = self.load_MC(base_path, recoil_species)[0]
        MC_drift = self.load_MC(base_path, recoil_species)[1] #only adding fiducial cuts
        tpcs = MC.keys()
        gain = {'iiwi': 1502, 'nene': 899, 'humu': 878, 'palila': 1033, 'tako': 807, 'elepaio': 797}
        W = 35.075
        tot_to_q = {}
        for tpc in tpcs:
            tot_to_q[tpc] = pd.DataFrame()
            tot_to_q[tpc]['tot_code'] = [i for i in range(0,14)]
        tot_to_q['iiwi']['conversion'] = [1833.00, 2345.17, 3017.33, 6001.54, 8891.71,
                                  11497.43, 14335.32, 18081.33, 22526.06, 27236.90,
                                  32056.16, 36955.09, 41874.75, 46794.40]

        tot_to_q['nene']['conversion'] = [1085.31, 2482.24, 4126.52, 5621.03, 7920.43,
                                  11667.35, 15117.97, 19489.23, 23211.63, 27483.98,
                                  32272.73, 37262.83, 42283.59, 47304.34]

        tot_to_q['humu']['conversion'] = [1758.29, 2324.41, 3679.37, 5433.43, 6862.72,
                                  10000.83, 13701.08, 17258.86, 21438.70, 25821.34,
                                  30153.82, 34460.74, 39042.80, 43624.85]

        tot_to_q['palila']['conversion'] = [1768.71, 2202.75, 2670.76, 4049.25, 6586.25,
                                    8954.45, 11551.60, 14428.46, 17618.81, 21140.34,
                                    24831.56, 28804.80, 33534.23, 40821.35]

        tot_to_q['tako']['conversion'] = [2761.20, 3077.66, 3509.80, 5475.02, 9230.59, 
                                  11955.00, 16837.46, 20761.78, 24514.73, 28445.96, 
                                  33071.27, 38033.29, 43011.21, 47989.15]

        tot_to_q['elepaio']['conversion'] = [1859.09, 2496.61, 4128.03, 6844.95, 9450.70,
                                     12158.68, 15125.31, 18507.89, 22166.14, 25826.40,
                                     29597.06, 33588.70, 38207.92, 42827.15]
        for tpc in tpcs:
            MC[tpc]['q_from_tot'] = MC[tpc]['tot']
            for i in range(0,len(MC[tpc])):
                MC[tpc]['q_from_tot'].iloc[i] = pd.Series(MC[tpc].iloc[i]['tot']).map(tot_to_q[tpc].set_index('tot_code')['conversion']).to_numpy()
            MC[tpc]['sumtot'] = [MC[tpc]['q_from_tot'][i].sum() for i in range(0,len(MC[tpc]))]
            MC[tpc]['reco_energy'] = MC[tpc]['sumtot']*35.075/gain[tpc]*1e-3
    
            #Add truth z fiducial cuts
            MC[tpc]['truth_z'] = [MC[tpc]['truth_center[2]'].iloc[i] for i in range(0,len(MC[tpc]))]
            MC[tpc] = MC[tpc].loc[(MC[tpc]['truth_z'] > zmin) & (MC[tpc]['truth_z'] < zmax)]
            MC[tpc].index = [i for i in range(0,len(MC[tpc]))]
            MC_drift[tpc]['truth_z'] = [MC_drift[tpc]['truth_center[2]'].iloc[i] for i in range(0,len(MC_drift[tpc]))]
            MC_drift[tpc] = MC_drift[tpc].loc[(MC_drift[tpc]['truth_z'] > zmin) & (MC_drift[tpc]['truth_z'] < zmax)]
            MC_drift[tpc].index = [i for i in range(0,len(MC_drift[tpc]))]


        return MC, MC_drift
            
    def add_saturation_fraction_and_mean_tot(self, base_path =  '~/data/phase3/simulation/resolution_paper/tpc_sims/', recoil_species = 'all', zmin = 0, zmax = 10):
        MC = self.apply_energy_calibrations(base_path = base_path, recoil_species = recoil_species, zmin = zmin, zmax = zmax)[0]
        tpcs = MC.keys()
        # mapping functions
        def get_saturation_fraction(dataframe):
            saturation_fraction = []
            for i in range(0,len(dataframe)):
                saturation_fraction.append(len([val for val in dataframe['tot'].iloc[i] if val == 13])/len(dataframe['tot'].iloc[i]))
            dataframe['saturation_fraction'] = saturation_fraction
            return dataframe
        def get_mean_tot(dataframe):
            tot_mean = []
            for i in range(0,len(dataframe)):
                tot_mean.append(dataframe['tot'].iloc[i].mean())
            dataframe['mean_tot'] = tot_mean
            return dataframe
        for tpc in tpcs:
            MC[tpc] = get_saturation_fraction(MC[tpc])
            MC[tpc] = get_mean_tot(MC[tpc])
        return MC

    def perform_saturation_corrections(self, base_path = '~/data/phase3/simulation/resolution_paper/tpc_sims/', recoil_species = 'all', zmin = 0, zmax = 10, poly_deg = 3):
        MC = self.add_saturation_fraction_and_mean_tot(base_path = base_path, recoil_species = recoil_species, zmin = zmin, zmax = zmax)
        tpcs = ['humu', 'iiwi', 'nene', 'palila', 'tako', 'elepaio']
        sat_frac_group = {} #dictionary of binned saturation fractions
        fit = {} #dictionary of polyfit calibration curve parameters for each TPC
        means = {} #dictionary of means in the saturation_frac group used to filter out NaNs
        sems = {} #dictionary of std errors in the saturation_frac group used to filter out NaNs
        for tpc in tpcs:
            #sat_frac_group[tpc] = MC[tpc].groupby(pd.cut(MC[tpc]['saturation_fraction'], bins = np.linspace(0.05,0.5,23)))
            sat_frac_group[tpc] = MC[tpc].groupby(pd.cut(MC[tpc]['saturation_fraction'], bins = np.linspace(0.0,0.45,46)))
            means[tpc] = sat_frac_group[tpc].mean().loc[sat_frac_group[tpc].sem()['reco_energy'].isna() == False]
            sems[tpc] = sat_frac_group[tpc].sem().loc[sat_frac_group[tpc].sem()['reco_energy'].isna() == False]
            #fit[tpc] = np.polyfit(means['humu']['saturation_fraction'], (means['humu']['reco_energy']/means['humu']['truth_energy']), poly_deg)
            fit[tpc] = np.polyfit(means[tpc]['saturation_fraction'], (means[tpc]['reco_energy']/means[tpc]['truth_energy']), poly_deg)
        #function for using fit calibration curve to correct for charge loss due to saturation
        def perform_saturation_correction(dataframe, fit):
            def get_correction_factor(x, fit):
                func = 0
                for i in range(0,len(fit)):
                    func += fit[i]*x**(len(fit)-i-1)
                return func
            dataframe['saturation_corrected_energy'] = 1/get_correction_factor(dataframe['saturation_fraction'], fit)*dataframe['reco_energy']
            return dataframe
        #apply function to MC
        for tpc in tpcs:
            MC[tpc] = perform_saturation_correction(MC[tpc], fit[tpc])
        return MC, fit #param [1] can be used with data

    def perform_threshold_corrections(self, base_path = '~/data/phase3/simulation/resolution_paper/tpc_sims/', recoil_species = 'all', zmin = 0, zmax = 10, poly_deg = 5):
        MC = self.perform_saturation_corrections(base_path, recoil_species, zmin, zmax, 3)[0]
        tpcs = MC.keys()
        mean_tot_group = {} #dictionary of binned mean_tots
        fit = {} #dictionary of polyfit calibration curve parameters for each TPC
        means = {} #dictionary of means in the saturation_frac group used to filter out NaNs
        sems = {} #dictionary of std errors in the saturation_frac group used to filter out NaNs
        for tpc in tpcs:
            mean_tot_group[tpc] = MC[tpc].groupby(pd.cut(MC[tpc]['mean_tot'], bins = np.linspace(0,8,26))) #cutoff at 8, may change later
            means[tpc] = mean_tot_group[tpc].mean().loc[mean_tot_group[tpc].sem()['saturation_corrected_energy'].isna() == False]
            sems[tpc] = mean_tot_group[tpc].sem().loc[mean_tot_group[tpc].sem()['saturation_corrected_energy'].isna() == False]
            fit[tpc] = np.polyfit(means[tpc]['mean_tot'], (means[tpc]['saturation_corrected_energy']/means[tpc]['truth_energy']), poly_deg)
        #function for applying threshold correction calibration curve
        def perform_threshold_correction(dataframe, fit):
            def get_correction_factor(x, fit):
                func = 0
                for i in range(0,len(fit)):
                    func += fit[i]*x**(len(fit)-i-1)
                return func
            dataframe['full_corrected_energy'] = 1/get_correction_factor(dataframe['mean_tot'], fit)*dataframe['saturation_corrected_energy']
            index = dataframe.loc[dataframe['mean_tot']>8].index.to_numpy()
            dataframe['full_corrected_energy'][index] = dataframe['saturation_corrected_energy'][index] #truncate full correction range
            return dataframe
        #apply function to MC
        for tpc in tpcs:
            MC[tpc] = perform_threshold_correction(MC[tpc], fit[tpc])
        return MC, fit #param[1] can be used with data

    def get_grouped_MC(self, base_path = '~/data/phase3/simulation/resolution_paper/tpc_sims/', recoil_species = 'all', zmin = 0, zmax = 10):
        after_threshold_data = self.apply_energy_calibrations(base_path = base_path, recoil_species = recoil_species, zmin = zmin, zmax = zmax)[0]
        IQF_data = self.apply_energy_calibrations(base_path = base_path, recoil_species = recoil_species, zmin = zmin, zmax = zmax)[1]
        tpcs = IQF_data.keys()
        IQF_group = {}
        after_threshold_group = {}
        for tpc in tpcs:
            IQF_group[tpc] = IQF_data[tpc].groupby(['truth_energy'])
            after_threshold_group[tpc] = after_threshold_data[tpc].groupby(['truth_energy'])
        return IQF_group, after_threshold_group

    def make_MC_summary_plots(self, base_path = '~/data/phase3/simulation/resolution_paper/tpc_sims/', recoil_species = 'all', zmin = 0, zmax = 10):
        MC, fit_thresh = self.perform_threshold_corrections(base_path, recoil_species, zmin, zmax, 5)
        #fit_sat2 = self.perform_saturation_corrections(base_path, recoil_species, zmin, zmax, 2)[1]
        fit_sat = self.perform_saturation_corrections(base_path, recoil_species, zmin, zmax, 3)[1]
        #fit_sat4 = self.perform_saturation_corrections(base_path, recoil_species, zmin, zmax, 4)[1]
        #fit_sat5 = self.perform_saturation_corrections(base_path, recoil_species, zmin, zmax, 5)[1]
        tpcs = MC.keys()
        
        def plot_fit(x, fit): #plots polyfit
            func = 0
            for i in range(0,len(fit)):
                func += fit[i]*x**(len(fit)-i-1)
            plt.plot(x,func, label='polyfit order %s'%(len(fit)-1))

        #Bin data for visualization of correction curves
        sat_frac_group = {}
        mean_tot_group = {}
        for tpc in tpcs:
            #sat_frac_group[tpc] = MC[tpc].groupby(pd.cut(MC[tpc]['saturation_fraction'], bins = np.linspace(0.05,0.5,23)))
            sat_frac_group[tpc] = MC[tpc].groupby(pd.cut(MC[tpc]['saturation_fraction'], bins = np.linspace(0.0,0.5,46)))
            mean_tot_group[tpc] = MC[tpc].groupby(pd.cut(MC[tpc]['mean_tot'], bins = np.linspace(0,8,26)))

        #Get grouped data for plotting ratio of ionization to truth energies vs truth energy
        after_thresh_group = {}
        IQF_group = self.get_grouped_MC(base_path = base_path, recoil_species = recoil_species, zmin = zmin, zmax = zmax)[0]
        for tpc in tpcs:
            after_thresh_group[tpc] = MC[tpc].groupby(['truth_energy'])
        #fig 1
        plt.figure(1, figsize = (10,10))
        i = 1
        x = np.linspace(0,1,101)
        for tpc in tpcs:
            plt.subplot(3,2,i)
            plt.plot(MC[tpc]['saturation_fraction'], MC[tpc]['reco_energy']/MC[tpc]['truth_energy'],'o', markersize = 1)
            #fit_df = MC[tpc].loc[MC[tpc]['saturation_fraction']>0.01]
            #fit2 = np.polyfit(fit_df['saturation_fraction'], fit_df['reco_energy']/fit_df['truth_energy'], 2)
            #fit3 = np.polyfit(fit_df['saturation_fraction'], fit_df['reco_energy']/fit_df['truth_energy'], 3)
            #fit4 = np.polyfit(fit_df['saturation_fraction'], fit_df['reco_energy']/fit_df['truth_energy'], 4)
            #fit5 = np.polyfit(fit_df['saturation_fraction'], fit_df['reco_energy']/fit_df['truth_energy'], 5)
            #fit6 = np.polyfit(fit_df['saturation_fraction'], fit_df['reco_energy']/fit_df['truth_energy'], 6)
            #plot_fit(x,fit2)
            #plot_fit(x,fit3)
            #plot_fit(x,fit4)
            #plot_fit(x,fit5)
            #plot_fit(x,fit6)
            plt.xlabel('Fraction of saturated pixels')
            plt.ylabel(r'$E_{reco}/E_{truth}$')
            plt.title(tpc)
            plt.xlim(0.0,1.1)
            plt.ylim(0,1)
            plt.grid()
            i += 1
        plt.tight_layout()
        plt.savefig('fig1.png')
        plt.clf()

        #fig 2
        means = {}
        sems = {}
        plt.figure(2, figsize = (10,10))
        i = 1
        gr = {}
        for tpc in tpcs:
            plt.subplot(3,2,i)
            means[tpc] = sat_frac_group[tpc].mean().loc[sat_frac_group[tpc].sem()['reco_energy'].isna() == False]
            sems[tpc] = sat_frac_group[tpc].sem().loc[sat_frac_group[tpc].sem()['reco_energy'].isna() == False]
            plt.errorbar(means[tpc]['saturation_fraction'], 
                 means[tpc]['reco_energy']/means[tpc]['truth_energy'],
                 sems[tpc]['reco_energy']/means[tpc]['truth_energy'],
                 sems[tpc]['saturation_fraction'], 'o')
            plot_fit(x,fit_sat[tpc])
            plt.xlabel('Fraction of saturated pixels')
            plt.ylabel(r'$E_{reco}/E_{truth}$')
            plt.title(tpc)
            plt.ylim(0,1)
            plt.xlim(-0.01,1.1)
            plt.grid()
            plt.legend()
            i += 1
        plt.tight_layout()
        plt.savefig('fig2.png')
        plt.clf()

        #fig 3
        plt.figure(3, figsize = (10,10))
        i = 1
        for tpc in tpcs:
            plt.subplot(3,2,i)
            plt.plot(MC[tpc]['mean_tot'], MC[tpc]['saturation_corrected_energy']/MC[tpc]['truth_energy'],'o',markersize = 1)
            plt.xlabel('Mean ToT per event')
            plt.ylabel(r'$E_{reco,sat cor}/E_{truth}$')
            plt.xlim(0,10)
            plt.ylim(0,1.5)
            plt.title(tpc)
            plt.grid()
            i += 1
        plt.tight_layout()
        plt.savefig('fig3.png')
        plt.clf()

        #fig 4
        plt.figure(4, figsize = (10,10))
        means = {}
        sems = {}
        i  = 1
        x = np.linspace(0,10,101)
        for tpc in tpcs:
            plt.subplot(3,2,i)
            means[tpc] = mean_tot_group[tpc].mean().loc[mean_tot_group[tpc].sem()['saturation_corrected_energy'].isna() == False]
            sems[tpc] = mean_tot_group[tpc].sem().loc[mean_tot_group[tpc].sem()['saturation_corrected_energy'].isna() == False]
            plt.errorbar(means[tpc]['mean_tot'], 
                 means[tpc]['saturation_corrected_energy']/means[tpc]['truth_energy'],
                 sems[tpc]['saturation_corrected_energy']/means[tpc]['truth_energy'],
                 sems[tpc]['mean_tot'], 'o')
            plot_fit(x, fit_thresh[tpc])
            plt.xlabel('Mean tot')
            plt.ylabel(r'$E_{reco}/E_{truth}$')
            plt.title(tpc)
            plt.ylim(-.5,1.5)
            plt.xlim(-0.01,8)
            plt.grid()
            i += 1
        plt.tight_layout()
        plt.savefig('fig4.png')
        plt.clf()

        #fig 5
        plt.figure(5, figsize = (10,10))
        i = 1
        for tpc in tpcs:
            plt.subplot(3,2,i)
            plt.plot(MC[tpc]['fit_length'], MC[tpc]['reco_energy'],'o',markersize = 1)
            plt.xlabel('Length [cm]')
            plt.ylabel(r'E_{reco}[keV]')
            
            #plt.ylim(0,1.5)
            plt.title(tpc)
            plt.grid()
            i += 1
        plt.tight_layout()
        plt.clf()

        #fig 6
        plt.figure(6, figsize = (10,10))
        i = 1
        for tpc in tpcs:
            plt.subplot(3,2,i)
            plt.plot(MC[tpc]['fit_length'], MC[tpc]['full_corrected_energy'],'o',markersize = 1)
            plt.xlabel('Length [cm]')
            plt.ylabel(r'E_{reco}[keV]')
            
            #plt.ylim(0,1.5)
            plt.title(tpc)
            plt.grid()
            i += 1
        plt.tight_layout()
        plt.clf()

        #fig 7
        gain = {'iiwi': 1502, 'nene': 899, 'humu': 878, 'palila': 1033, 'tako': 807, 'elepaio': 797}
        W = 35.075
        
        plt.rc('legend', fontsize=8)
        plt.rc('xtick', labelsize=12)
        plt.rc('ytick', labelsize=12)
        plt.rc('axes', labelsize=14)
        plt.rc('axes', titlesize=14)
        plt.figure(figsize = (12,12))
        i = 1
        for tpc in tpcs:
            plt.subplot(3,2,i)
            plt.errorbar(IQF_group[tpc].mean().index,
                 IQF_group[tpc].mean()['truth_charge']*W*1e-03/IQF_group[tpc].mean().index,
                 IQF_group[tpc].sem()['truth_charge']*W*1e-03/IQF_group[tpc].mean().index,
                 [0 for i in range(0,len(IQF_group[tpc].mean()))],
                 'o',markersize = 3, label = r'After quenching', color = 'red')

            plt.errorbar(after_thresh_group[tpc].mean().index,
                 (after_thresh_group[tpc].mean()['qsum']*W*1e-03/gain[tpc])/after_thresh_group[tpc].mean().index,
                 (after_thresh_group[tpc].sem()['qsum']*W*1e-03/gain[tpc])/after_thresh_group[tpc].mean().index,
                 [0 for i in range(0,len(after_thresh_group[tpc].mean()))],
                 'o',markersize = 3, label = r'After threshold loss', color = 'magenta')

            plt.errorbar(after_thresh_group[tpc].mean().index,
                 (after_thresh_group[tpc].mean()['sumtot']*W*1e-03/gain[tpc])/after_thresh_group[tpc].mean().index,
                 (after_thresh_group[tpc].sem()['sumtot']*W*1e-03/gain[tpc])/after_thresh_group[tpc].mean().index,
                 [0 for i in range(0,len(after_thresh_group[tpc].mean()))],
                 'o',markersize = 3, label = r'After saturation loss', color = 'blue')

            
            plt.errorbar(after_thresh_group[tpc].mean().index,
                 (after_thresh_group[tpc].mean()['saturation_corrected_energy'])/after_thresh_group[tpc].mean().index,
                 (after_thresh_group[tpc].sem()['saturation_corrected_energy'])/after_thresh_group[tpc].mean().index,
                 [0 for i in range(0,len(after_thresh_group[tpc].mean()))],
                 'o',markersize = 3, label = r'Saturation correction', color = 'indigo')

            plt.errorbar(after_thresh_group[tpc].mean().index,
                 (after_thresh_group[tpc].mean()['full_corrected_energy'])/after_thresh_group[tpc].mean().index,
                 (after_thresh_group[tpc].sem()['full_corrected_energy'])/after_thresh_group[tpc].mean().index,
                 [0 for i in range(0,len(after_thresh_group[tpc].mean()))],
                 'o',markersize = 3, label = r'Threshold and saturation correction', color = 'green')

            plt.ylabel(r'$E_{ion}/E_{truth}$')
            plt.xlabel(r'$E_{truth}$ [keV]')
            plt.ylim(0,1.2)
            plt.xlim(0,1050)
            plt.legend(loc ='lower right')
            plt.grid()
            plt.title(tpc)
            i += 1
        plt.tight_layout()
        plt.savefig("/home/jeef/Pictures/truth_frac_%s_zmin-%s_zmax_%s.png"%(recoil_species, zmin, zmax))
        plt.clf()

        # fig 8: Resolutions
        h_thresh_sat = {}
        h_thresh = {}
        h_quench = {}
        h_full_cor = {}
        h_sat_cor = {}
        vals_thresh_sat = {}
        vals_thresh = {}
        vals_quench = {}
        vals_full_cor = {}
        vals_sat_cor = {}
        rms_errors_quench = {}
        rms_errors_thresh = {}
        rms_errors_thresh_sat = {}
        rms_errors_full_cor = {}
        rms_errors_sat_cor = {}
        reso_errors_thresh_sat = {}
        reso_errors_thresh = {}
        reso_errors_quench = {}
        reso_errors_full_cor = {}
        reso_errors_sat_cor = {}
        
        IQF_data = self.apply_energy_calibrations(base_path = base_path, recoil_species = recoil_species, zmin = zmin, zmax = zmax)[1]

        for tpc in tpcs:
            rms_errors_sat_cor[tpc] = []
            for energy in MC[tpc]['truth_energy'].unique(): #sat_cor loop
                h_sat_cor[tpc+'_'+str(energy)] = ROOT.TH1F('%s_%s'%(tpc,energy), '%s_%s'%(tpc,energy), 21, 0, 100)
                vals_sat_cor[tpc+'_'+str(energy)] = array.array('d', after_thresh_group[tpc].mean().loc[MC[tpc]['truth_energy'] == energy]['saturation_corrected_energy'])
                for i in range(0,len(vals_sat_cor[tpc+'_'+str(energy)])):
                    h_sat_cor[tpc+'_'+str(energy)].Fill(vals_sat_cor[tpc+'_'+str(energy)][i])
                rms_errors_sat_cor[tpc].append(h_sat_cor[tpc+'_'+str(energy)].GetRMSError())
            rms_errors_sat_cor[tpc] = np.array(rms_errors_sat_cor[tpc])
            reso_errors_sat_cor[tpc] = after_thresh_group[tpc].std()['saturation_corrected_energy']/after_thresh_group[tpc].mean()['saturation_corrected_energy'] * np.sqrt((rms_errors_sat_cor[tpc]/(after_thresh_group[tpc].std()['saturation_corrected_energy']))**2 + (after_thresh_group[tpc].std()['saturation_corrected_energy']/(after_thresh_group[tpc].mean()['saturation_corrected_energy']))**2)
            
        for tpc in tpcs:
            rms_errors_full_cor[tpc] = []
            for energy in MC[tpc]['truth_energy'].unique(): #full_cor loop
                h_full_cor[tpc+'_'+str(energy)] = ROOT.TH1F('%s_%s'%(tpc,energy), '%s_%s'%(tpc,energy), 21, 0, 100)
                vals_full_cor[tpc+'_'+str(energy)] = array.array('d', after_thresh_group[tpc].mean().loc[MC[tpc]['truth_energy'] == energy]['full_corrected_energy'])
                for i in range(0,len(vals_full_cor[tpc+'_'+str(energy)])):
                    h_full_cor[tpc+'_'+str(energy)].Fill(vals_full_cor[tpc+'_'+str(energy)][i])
                rms_errors_full_cor[tpc].append(h_full_cor[tpc+'_'+str(energy)].GetRMSError())
            rms_errors_full_cor[tpc] = np.array(rms_errors_full_cor[tpc])
            reso_errors_full_cor[tpc] = after_thresh_group[tpc].std()['full_corrected_energy']/after_thresh_group[tpc].mean()['full_corrected_energy'] * np.sqrt((rms_errors_sat_cor[tpc]/(after_thresh_group[tpc].std()['full_corrected_energy']))**2 + (after_thresh_group[tpc].std()['full_corrected_energy']/(after_thresh_group[tpc].mean()['full_corrected_energy']))**2)
            
        for tpc in tpcs:
            rms_errors_quench[tpc] = []
            rms_errors_thresh[tpc] = []
            rms_errors_thresh_sat[tpc] = []
            for energy in IQF_data[tpc]['truth_energy'].unique(): #quenching loop
                h_quench[tpc+'_'+str(energy)] = ROOT.TH1F('%s_%s'%(tpc,energy), '%s_%s'%(tpc,energy), 21, 0, 100)
                vals_quench[tpc+'_'+str(energy)] = array.array('d', after_thresh_group[tpc].mean().loc[MC[tpc]['truth_energy'] == energy]['truth_charge']*W*1e-03)
                for i in range(0,len(vals_quench[tpc+'_'+str(energy)])):
                    h_quench[tpc+'_'+str(energy)].Fill(vals_quench[tpc+'_'+str(energy)][i])
                rms_errors_quench[tpc].append(h_quench[tpc+'_'+str(energy)].GetRMSError())
            for energy in MC[tpc]['truth_energy'].unique(): #thresh and saturation
                h_thresh_sat[tpc+'_'+str(energy)] = ROOT.TH1F('%s_%s'%(tpc,energy), '%s_%s'%(tpc,energy), 21, 0, 100)
                h_thresh[tpc+'_'+str(energy)] = ROOT.TH1F('%s_%s'%(tpc,energy), '%s_%s'%(tpc,energy), 21, 0, 100)
                vals_thresh_sat[tpc+'_'+str(energy)] = array.array('d', after_thresh_group[tpc].mean().loc[MC[tpc]['truth_energy'] == energy]['sumtot']*W*1e-03/gain[tpc])
                vals_thresh[tpc+'_'+str(energy)] = array.array('d', after_thresh_group[tpc].mean().loc[MC[tpc]['truth_energy'] == energy]['qsum']*W*1e-03/gain[tpc])
                for i in range(0,len(vals_thresh_sat[tpc+'_'+str(energy)])):
                    h_thresh_sat[tpc+'_'+str(energy)].Fill(vals_thresh_sat[tpc+'_'+str(energy)][i])
                for i in range(0,len(vals_thresh[tpc+'_'+str(energy)])):
                    h_thresh[tpc+'_'+str(energy)].Fill(vals_thresh[tpc+'_'+str(energy)][i])
                rms_errors_thresh_sat[tpc].append(h_thresh_sat[tpc+'_'+str(energy)].GetRMSError())
                rms_errors_thresh[tpc].append(h_thresh[tpc+'_'+str(energy)].GetRMSError())
            rms_errors_thresh_sat[tpc] = np.array(rms_errors_thresh_sat[tpc])
            rms_errors_thresh[tpc] = np.array(rms_errors_thresh[tpc])
            rms_errors_quench[tpc] = np.array(rms_errors_quench[tpc])
            reso_errors_quench[tpc] = IQF_group[tpc].std()['truth_charge']*W*1e-03/IQF_group[tpc].mean()['truth_charge']*W*1e-03 * np.sqrt((rms_errors_quench[tpc]/IQF_group[tpc].std()['truth_charge']*W*1e-03)**2 + (IQF_group[tpc].std()['truth_charge']*W*1e-03/IQF_group[tpc].mean()['truth_charge']*W*1e-03)**2)
            reso_errors_thresh[tpc] = (after_thresh_group[tpc].std()['qsum']*W*1e-03/gain[tpc])/(after_thresh_group[tpc].mean()['qsum']*W*1e-03/gain[tpc]) * np.sqrt((rms_errors_thresh[tpc]/(after_thresh_group[tpc].std()['qsum']*W*1e-03)/gain[tpc])**2 + ((after_thresh_group[tpc].std()['qsum']*W*1e-03/gain[tpc])/(after_thresh_group[tpc].mean()['qsum']*W*1e-03/gain[tpc]))**2)
            reso_errors_thresh_sat[tpc] = (after_thresh_group[tpc].std()['sumtot']*W*1e-03/gain[tpc])/(after_thresh_group[tpc].mean()['sumtot']*W*1e-03/gain[tpc]) * np.sqrt((rms_errors_thresh_sat[tpc]/(after_thresh_group[tpc].std()['sumtot']*W*1e-03)/gain[tpc])**2 + ((after_thresh_group[tpc].std()['sumtot']*W*1e-03/gain[tpc])/(after_thresh_group[tpc].mean()['sumtot']*W*1e-03/gain[tpc]))**2)

        # Make plots
        plt.rc('legend', fontsize=8)
        plt.rc('xtick', labelsize=12)
        plt.rc('ytick', labelsize=12)
        plt.rc('axes', labelsize=14)
        plt.rc('axes', titlesize=14)
        plt.figure(figsize = (12,12))

        i = 1
        for tpc in tpcs:
            plt.subplot(3,2,i)
            plt.errorbar(IQF_group[tpc].mean().index, IQF_group[tpc].std()['truth_charge']*35.075*1e-03/IQF_group[tpc].mean()['truth_charge']*35.075*1e-03, reso_errors_quench[tpc], [0 for i in range(0,len(reso_errors_quench[tpc]))], 'o',markersize = 3, label = r'After quenching', color = 'red', lw = 1)
            plt.errorbar(after_thresh_group[tpc].mean().index, (after_thresh_group[tpc].std()['qsum']*35.075*1e-03/gain[tpc])/(after_thresh_group[tpc].mean()['qsum']*35.075*1e-03/gain[tpc]), reso_errors_thresh[tpc], [0 for i in range(0,len(reso_errors_thresh[tpc]))], 'o',markersize = 3, label = r'After threshold loss', color = 'magenta', lw = 1)
            plt.errorbar(after_thresh_group[tpc].mean().index, (after_thresh_group[tpc].std()['sumtot']*35.075*1e-03/gain[tpc])/(after_thresh_group[tpc].mean()['sumtot']*35.075*1e-03/gain[tpc]), reso_errors_thresh_sat[tpc], [0 for i in range(0,len(reso_errors_thresh_sat[tpc]))], 'o',markersize = 3, label = r'After saturation loss', color = 'blue', lw = 1)
            plt.errorbar(after_thresh_group[tpc].mean().index, (after_thresh_group[tpc].std()['saturation_corrected_energy'])/(after_thresh_group[tpc].mean()['saturation_corrected_energy']), reso_errors_sat_cor[tpc], [0 for i in range(0,len(reso_errors_sat_cor[tpc]))], 's',markersize = 3, label = r'Sat. loss correction', color = 'indigo', lw = 1)
            plt.errorbar(after_thresh_group[tpc].mean().index, (after_thresh_group[tpc].std()['full_corrected_energy'])/(after_thresh_group[tpc].mean()['full_corrected_energy']), reso_errors_full_cor[tpc], [0 for i in range(0,len(reso_errors_full_cor[tpc]))], 's',markersize = 3, label = r'Sat + thresh correction', color = 'green', lw = 1)
            plt.xlabel(r'$E_{truth}$ [keV]')
            plt.ylabel(r'$\sigma_E/E$')
            plt.xlim(0,1020)
            plt.ylim(0,0.5)
            plt.title(tpc)
            plt.grid()
            plt.legend()
            i += 1
        plt.tight_layout()
        plt.savefig("/home/jeef/Pictures/resolution_%s_zmin-%s_zmax_%s.png"%(recoil_species, zmin, zmax))
        plt.clf()
        
            
class tpc_calibration(simulation):

    def __init__(self, corrected_energy = False):
        #super().__init__() # inherit simulation classes methods
        self.alphas, self.scale_factor, self.scale_factor_err = self.calibrate_alphas(corrected_energy = corrected_energy)
        self.recoils = self.calibrate_recoils(corrected_energy = corrected_energy)
        self.updated_recoils = self.determine_new_angles_and_lengths()
        self.write_to_root_file()
        #pass
        
    def get_tpc_list(self, tpc_list = ['iiwi', 'humu', 'nene', 'tako', 'palila', 'elepaio']):
        return tpc_list

    def load_alphas(self):
        tpcs = self.get_tpc_list()
        df = {}
        for tpc in tpcs:
            df[tpc] = ur.open('~/data/phase3/spring_2020/maintenance_day_test/%s_alphas_all.root'%(tpc))[ur.open('~/data/phase3/spring_2020/maintenance_day_test/%s_alphas_all.root'%(tpc)).keys()[0]].pandas.df(flatten=False)
            df[tpc]['BCID_range'] = [df[tpc]['BCID'][i].max()-df[tpc]['BCID'][i].min() for i in range (0,len(df[tpc]))]
            df[tpc] = df[tpc].loc[(df[tpc]['track_energy'] > 400) & (df[tpc]['BCID_range'] < 90)] #selection for alphas
            df[tpc] = df[tpc].loc[(np.abs(df[tpc]['phi']) < 5) | (np.abs(df[tpc]['phi']-180) < 5) | (np.abs(-180 - df[tpc]['phi']) < 5)]
            df[tpc] = df[tpc].loc[(df[tpc]['theta']>85) & (df[tpc]['theta']<95)]
            df[tpc]['track_energy'] = df[tpc]['track_charge']*35.075/2000*1e-3 # "uncalibrates" energy to start fresh
            df[tpc]['pixel_energy'] = df[tpc]['pixel_charge']*35.075/2000*1e-3 # "uncalibrates" energy to start fresh
        return df #returns dictionary of horizontal alphas

    def get_study_data(self):
        tpcs = self.get_tpc_list()
        df = {}
        for tpc in tpcs:
            df[tpc] = ur.open('~/data/phase3/spring_2020/05-09-20/tpc_root_files/%s_all_new.root'%(tpc))[ur.open('~/data/phase3/spring_2020/05-09-20/tpc_root_files/%s_all_new.root'%(tpc)).keys()[0]].pandas.df(flatten=False)
            df[tpc]['track_energy'] = df[tpc]['track_charge']*35.075/2000*1e-3 # "uncalibrates" energy to start fresh
            df[tpc]['pixel_energy'] = df[tpc]['pixel_charge']*35.075/2000*1e-3 # "uncalibrates" energy to start fresh
        return df

    def loose_neutron_cuts(self, dataframe): #computes minimal set of cuts to start
        dataframe = dataframe.loc[(dataframe['hitside_row_min']==0) & (dataframe['hitside_row_max']==0) & (dataframe['hitside_col_min']==0) & (dataframe['hitside_col_max']==0)]
        return dataframe

    def load_recoils(self):
        df = self.get_study_data()
        tpcs = df.keys()
        for tpc in tpcs:
            df[tpc] = self.loose_neutron_cuts(df[tpc])
        return df

    #def get_saturation_correction_factor(self, x, tpc, fit, recoil_species = 'He'):
    #    #return -1.03*x**4 +3.73*x**3 -4.65*x**2 + 0.94*x + 0.79 #Michael's correction factor
    # Perform fit on MC using simulation class
    #    func = 0
    #    for i in range(0,len(fit)):
    #        func += fit[i]*x**(len(fit)-i-1)
    #    return func
        
    def get_correction_factor(self, x, tpc, fit):
        #return 6.7e-5*x**5 -7.9e-4*x**4 - 3.5e-3 * x**3 + 4.9e-2*x**2 +7e-2*x + 1.2e-1 #Michael
        func = 0
        for i in range(0,len(fit)):
            func += fit[i]*x**(len(fit)-i-1)
        return func
        
    def perform_saturation_correction(self, dataframe, tpc, fit):
        saturation_fraction = []
        for i in range(0,len(dataframe)):
            saturation_fraction.append(len([val for val in dataframe['tot'].iloc[i] if val == 13])/len(dataframe['tot'].iloc[i]))
        dataframe['saturation_fraction'] = saturation_fraction
        dataframe['saturation_corrected_energy'] = 1/self.get_correction_factor(dataframe['saturation_fraction'], tpc, fit)*dataframe['track_energy']
        return dataframe

    def perform_threshold_correction(self, dataframe, tpc, fit):
        tot_mean = []
        for i in range(0,len(dataframe)):
            tot_mean.append(dataframe['tot'].iloc[i].mean())
        dataframe['mean_tot'] = tot_mean
        dataframe['full_corrected_energy'] = 1/self.get_correction_factor(dataframe['mean_tot'], tpc, fit)*dataframe['saturation_corrected_energy']
        index = dataframe.loc[dataframe['mean_tot']>8].index.to_numpy()
        dataframe['full_corrected_energy'][index] = dataframe['saturation_corrected_energy'][index] #truncate full correction range
        return dataframe

    def correct_alphas(self): #corrects alphas for saturation and threshold
        alphas = self.load_alphas()
        tpcs = alphas.keys()
        fit_sat = self.perform_saturation_corrections(recoil_species = 'He', poly_deg = 3)[1] #dict w/ tpc names as keys. Comes from simulation class
        fit_thresh = self.perform_threshold_corrections(recoil_species = 'He', poly_deg = 5)[1] #dict w/ tpc names as keys. From simulation class
        for tpc in tpcs:
            alphas[tpc] = self.perform_saturation_correction(alphas[tpc], tpc, fit_sat[tpc])
            alphas[tpc] = self.perform_threshold_correction(alphas[tpc], tpc, fit_thresh[tpc])
        return alphas

    def correct_recoils(self): #corrects alphas for saturation and threshold
        recoils = self.load_recoils()
        tpcs = recoils.keys()
        fit_sat = self.perform_saturation_corrections(recoil_species = 'He', poly_deg = 3)[1] #dict w/ tpc names as keys
        fit_thresh = self.perform_threshold_corrections(recoil_species = 'He', poly_deg = 5)[1] #dict w/ tpc names as keys
        for tpc in tpcs:
            recoils[tpc] = self.perform_saturation_correction(recoils[tpc], tpc, fit_sat[tpc])
            recoils[tpc] = self.perform_threshold_correction(recoils[tpc], tpc, fit_thresh[tpc])
        return recoils

    def calibrate_alphas(self, corrected_energy): #calibrates alphas to a de/dx reference value of 500 keV/cm
        if corrected_energy == True:
            alphas = self.correct_alphas()
            ekey = 'full_corrected_energy'
        else:
            alphas = self.load_alphas()
            ekey = 'track_energy'
        tpcs = alphas.keys()
        scale_factor = {}
        scale_factor_err = {}
        for tpc in tpcs:
            #alphas[tpc]['dedx'] = alphas[tpc][ekey]/alphas[tpc][lkey]
            scale_factor[tpc] = 1430/((alphas[tpc][ekey]*np.sin(alphas[tpc]['theta']*np.pi/180)).mean())
            scale_factor_err[tpc] = scale_factor[tpc]*(alphas[tpc][ekey]*np.sin(alphas[tpc]['theta']*np.pi/180)).sem()/((alphas[tpc][ekey]*np.sin(alphas[tpc]['theta']*np.pi/180)).mean())
            #print((alphas[tpc][ekey]*np.sin(alphas[tpc]['theta']*np.pi/180)).sem())
            alphas[tpc][ekey] = scale_factor[tpc]*alphas[tpc][ekey]
            alphas[tpc][ekey+'_err'] = scale_factor[tpc]*alphas[tpc][ekey]*(alphas[tpc][ekey]*np.sin(alphas[tpc]['theta']*np.pi/180)).sem()/((alphas[tpc][ekey]*np.sin(alphas[tpc]['theta']*np.pi/180)).mean())
            #print(tpc, scale_factor_err[tpc]/scale_factor[tpc])
        return alphas, scale_factor, scale_factor_err

    def calibrate_recoils(self, corrected_energy): #calibrates recoils to a de/dx reference value of 500 keV/cm
        if corrected_energy == True:
            recoils = self.correct_recoils()
            ekey = 'full_corrected_energy'
        else:
            recoils = self.load_recoils()
            ekey = 'track_energy'
        scale_factor, scale_factor_err = self.scale_factor, self.scale_factor_err
        tpcs = recoils.keys()
        for tpc in tpcs:
            recoils[tpc]['track_energy'] = scale_factor[tpc]*recoils[tpc]['track_energy']
            recoils[tpc]['track_energy_err'] = scale_factor_err[tpc]*recoils[tpc]['track_energy']
            #recoils[tpc]['track_energy_err'] = scale_factor_err[tpc]*alphas[tpc]['track_energy'].mean()
            #recoils[tpc]['track_energy_err'] = scale_factor_err[tpc]*recoils[tpc]['track_energy']/scale_factor[tpc]
            #print(tpc, recoils[tpc]['track_energy_err']/recoils[tpc]['track_energy'])
            if corrected_energy == True:
                recoils[tpc]['saturation_corrected_energy'] = scale_factor[tpc]*recoils[tpc]['saturation_corrected_energy']
                recoils[tpc]['full_corrected_energy'] = scale_factor[tpc]*recoils[tpc]['full_corrected_energy']
        return recoils
            
    def load_colors(self): #makes colorblind friendly palette
        color_codes = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00'] #friendly colors for color blind
        #color_codes = ['#FFD700', '#4B0082', '#FF00FF', '#00FFFF', '#FF0000', '#808080']
        return color_codes

    def plot_alpha_distributions(self, corrected_energy = True, calibrated = True, gain = False):
        if calibrated == True:
            alphas, scale_factor, scale_factor_err = self.calibrate_alphas(corrected_energy)
        else:
            alphas = self.correct_alphas()
        if corrected_energy == True:
            ekey = 'full_corrected_energy'
        else:
            ekey = 'track_energy'
        tpcs = ['elepaio', 'tako', 'palila', 'iiwi', 'nene', 'humu']
        locations = ['z=-14m', 'z=-8.0m', 'z=-5.6m', 'z=+6.6m', 'z=+14m', 'z=+16m']
        color_codes = self.load_colors()
        #colors = [color_codes[i] for i in range(0,len(tpcs))]
        colors = ['tab:blue' for i in range(0,len(tpcs))]
        i = 0
        for tpc in tpcs:
            #plt.hist(alphas[tpc][ekey]/alphas[tpc][lkey]*1e4, bins = 30, range = (200, 800), linewidth = 2, edgecolor = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=1), color = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=0.2), label = tpc)
            if gain == True:
                alphas[tpc][ekey] = alphas[tpc]['track_charge']*35.075/alphas[tpc][ekey]*1e-3
                #print(tpc, alphas[tpc][ekey].unique(), alphas[tpc][ekey+'_err'].unique())
                #print(tpc, alphas[tpc][ekey].mean(), scale_factor_err[tpc]*alphas[tpc][ekey].mean())
                plt.bar(locations[i], alphas[tpc][ekey].mean(), yerr = scale_factor_err[tpc]*alphas[tpc][ekey].mean(), linewidth = 2, edgecolor = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=1), color = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=0.2))
            elif calibrated == True:
                plt.hist(alphas[tpc][ekey], bins = 60, range = (0, 2000), linewidth = 2, edgecolor = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=1), color = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=0.2), label = '%s'%(tpc)) #hatch = '/\\'
            elif calibrated == False:
                plt.hist(alphas[tpc][ekey], bins = 60, range = (0, 2000), linewidth = 2, edgecolor = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=1), color = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=0.2), label = '%s uncalibrated'%(tpc))
            i += 1
        if gain == False:
            plt.xlabel(r'$E\cos(\theta)$ [keV]') #theta is z-x angle (equiv to sine of theta in our dataframes)
        #plt.title(title)
        #plt.savefig(fig)
        #plt.clf()
        plt.show()
    def plot_calibrated_recoil_distributions(self):
        recoils = self.recoils
        if corrected_energy == True:
            ekey = 'full_corrected_energy'
        else:
            ekey = 'track_energy'
        tpcs = recoils.keys()
        color_codes = self.load_colors()
        colors = [color_codes[i] for i in range(0,len(tpcs))]
        i = 0
        for tpc in tpcs:
            plt.plot(recoils[tpc][lkey], recoils[tpc][ekey], 'o', markersize = 1, alpha = 0.5, label = tpc)
            plt.legend()
            plt.ylabel('Energy [keV]')
            plt.xlabel('length [um]')
            plt.xlim(0,10000)
            plt.ylim(0,200)
            #plt.savefig('_%s.png'%(i))
            i += 1
        #plt.clf()
        plt.show()

        
    def determine_new_angles_and_lengths(self):
        recoils = self.recoils
        tpcs = recoils.keys()
        ### Scale z by correct drift speed ###
        for tpc in tpcs:
            if tpc != 'elepaio':
                recoils[tpc]['z'] = recoils[tpc]['z']*7/8
            else:
                recoils[tpc]['z'] = recoils[tpc]['z']*5.5/8
        ### PERFORM SVD FITTER ###
        def fitsvd_numpy(df, i): #faster than root fit, so this is standard
            x = df.iloc[i]['x']
            y = df.iloc[i]['y']
            z = df.iloc[i]['z']
            q = df.iloc[i]['pixel_charge']
            data = np.concatenate((x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis]), axis=1)
            datamean = data.mean(axis=0)
            uu, dd, vv = np.linalg.svd(data - datamean)
            projection = []
            for point in data:
                projection += [np.dot(point, vv[0])]
            vec = TVector3(vv[0][0],vv[0][1],vv[0][2])
            maxp = max(projection)
            minp = min(projection)
            length = maxp - minp    
            midp = 0.5*float(maxp+minp)
            head_charge = 0
            tail_charge = 0
            i=0
            for p in projection:
                if p > midp:
                    head_charge += q[i]
                else:
                    tail_charge += q[i]
                i += 1
            head_charge_fraction = head_charge/(head_charge+tail_charge)
            theta = vec.Theta() #RADIANS
            phi = vec.Phi() #RADIANS

            if phi < -np.pi/2 or phi > np.pi/2: #fold phi to be restricted from 0 to pi
                vec = -1*vec
                theta = vec.Theta()
                phi = vec.Phi()

            return length, theta, phi, head_charge_fraction, head_charge, tail_charge

        def update_dataframe(df):
            length = []
            theta = []
            phi = []
            head_charge_frac = []
            head_charge = []
            tail_charge = []
            for j in range(0,len(df)):
                l, t, p, q, hc, tc = fitsvd_numpy(df, j)
                length.append(l)
                theta.append(t)
                phi.append(p)
                head_charge_frac.append(q)
                head_charge.append(hc)
                tail_charge.append(tc)
            df['new_length'] = length
            df['new_theta'] = theta
            df['new_phi'] = phi
            df['new_head_q_frac'] = head_charge_frac
            df['new_head_q'] = head_charge
            df['new_tail_q'] = tail_charge
            return df

        for tpc in tpcs:
            recoils[tpc] = update_dataframe(recoils[tpc])
            
        return recoils
        
    def apply_recoil_cuts(self): #Cuts without xray veto
        recoils = self.updated_recoils
        x = {'iiwi': np.array([850, 2500, 15000]), 'nene': np.array([950, 2500, 16000]), 'humu': np.array([1950, 3000, 20000]),
             'palila': np.array([900, 2350, 15000]), 'tako': np.array([1400, 2500, 15000]), 'elepaio': np.array([1050, 2500, 15000])}
        y = np.array([6,20,800])
        cut = {}
        tpcs = recoils.keys()
        for tpc in tpcs:
            cut[tpc] = np.polyfit(x[tpc],y,2)
            recoils[tpc]['is_recoil'] = 0
            index = recoils[tpc].loc[(recoils[tpc]['track_energy']>=(cut[tpc][0]*recoils[tpc]['new_length']**2 + cut[tpc][1]*recoils[tpc]['new_length'] + cut[tpc][2]))].index.to_numpy()
            recoils[tpc]['is_recoil'][index]=1
        
        ''' ###old###
        recoils =  self.recoils
        y = np.array([6,20,800])
        for tpc in recoils.keys():
            #recoils[tpc] = recoils[tpc].loc[recoils[tpc]['track_energy']>8]
            if tpc == 'iiwi':
                x = np.array([1200, 1900, 15000])
            elif tpc == 'humu':
                x = np.array([1950, 3000, 20000])
            elif tpc == 'nene':
                x = np.array([950, 1900, 15000])
            elif tpc == 'tako':
                x = np.array([1000, 1900, 15000])
            elif tpc == 'palila':
                x = np.array([1000, 1750, 15000])
            else:
                x = np.array([1050, 2000, 15000])
            cut = np.polyfit(x,y,2)
            recoils[tpc] = recoils[tpc].loc[recoils[tpc]['track_energy'] > cut[0]*recoils[tpc]['length']**2+cut[1]*recoils[tpc]['length']+cut[2]] #after this cut, only recoil bands remain
        '''
        return recoils

    def write_to_root_file(self, outdir = '~/data/phase3/spring_2020/05-09-20/tpc_root_files/'):
        #if recoils_only == False:
        #    recoils = self.recoils
        #else:
        recoils = self.apply_recoil_cuts()
        tpcs = recoils.keys()
        #tpcs = ['nene','humu','tako','palila','elepaio']
        for tpc in tpcs:            
            recoils[tpc].index = [i for i in range(0,len(recoils[tpc]))]
            recoils[tpc]['event_number'] = recoils[tpc].index
            keys = [val for val in recoils[tpc].columns]
            #if recoils_only == True:
            #    output = ROOT.TFile(outdir + '%s_all_recoils_only_even_newester2.root'%(tpc), 'new')
            #else:
            output = ROOT.TFile(outdir + '%s_all_newester4.root'%(tpc), 'recreate')
            tout = ROOT.TTree('data','data')
            branches = {}
            data={}
        
            for key in keys:
                if recoils[tpc][key].dtype == "O": #Determines the size of an array in a dataframe to be pushed to ntuple
                    data[key]=array.array('d',[0 for j in range(0,50000)])
                    branches[key]=tout.Branch("%s"%(key), data[key], "%s[npoints]/D"%(key))
                elif key == 'npoints':
                    data[key]=array.array('i',[0])
                    branches[key]=tout.Branch("%s"%(key), data[key], "%s/I"%(key))
                else:
                    data[key]=array.array('d',[0])
                    branches[key]=tout.Branch("%s"%(key), data[key], "%s/D"%(key))

            for j in range(0,len(recoils[tpc])):
                data['npoints'][0] = recoils[tpc]['npoints'].to_numpy()[j].astype(int)
                for key in keys:
                    if recoils[tpc][key].dtype == "O":
                        for i in range(0,data['npoints'][0]):
                            data[key][i]=recoils[tpc][key][j][i]
                    elif key != 'npoints':
                        data[key][0]=recoils[tpc][key][j]
                print("Filling event %s out of %s"%(j+1, len(recoils[tpc])))
                tout.Fill()

            output.Write()
            output.Close()
            
         
t = tpc_calibration()

#Old code below#

'''
def plot_uncalibrated_recoil_distributions(self):
        recoils = self.load_recoils()
        tpcs = recoils.keys()
        color_codes = self.load_colors()
        colors = [color_codes[i] for i in range(0,len(tpcs))]
        i = 0
        for tpc in tpcs:
            plt.plot(recoils[tpc]['length'], recoils[tpc]['track_energy'], 'o', markersize = 1, alpha = 0.5, label = tpc)
            plt.legend()
            plt.ylabel('Energy [keV]')
            plt.xlabel('length [um]')
            plt.xlim(0,20000)
            plt.ylim(0,1000)
            plt.savefig('calibration_figures/uncalibrated_EvL_%s.png'%(i))
            i += 1
        plt.clf()

    def plot_corrected_uncalibrated_alpha_distributions(self):
        alphas = self.correct_alphas()
        tpcs = alphas.keys()
        color_codes = self.load_colors()
        colors = [color_codes[i] for i in range(0,len(tpcs))]
        i = 0
        for tpc in tpcs:
            plt.hist(alphas[tpc]['full_corrected_energy']/alphas[tpc]['length']*1e4, bins = 30, range = (200, 800), linewidth = 2, edgecolor = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=1), color = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=0.2), label = tpc)
            i += 1
        plt.xlabel('dE/dx [keV/cm]')
        plt.legend()
        plt.title('Corrected Uncalibrated Horizontal Alpha Distribution')
        plt.savefig('calibration_figures/corrected_uncalibrated_alphas.png')
        plt.clf()

    def plot_corrected_uncalibrated_recoil_distributions(self):
        recoils = self.correct_recoils()
        tpcs = recoils.keys()
        color_codes = self.load_colors()
        colors = [color_codes[i] for i in range(0,len(tpcs))]
        i = 0
        for tpc in tpcs:
            plt.plot(recoils[tpc]['length'], recoils[tpc]['full_corrected_energy'], 'o', markersize = 1, alpha = 0.5, label = tpc)
            plt.legend()
            plt.ylabel('Energy [keV]')
            plt.xlabel('length [um]')
            plt.xlim(0,20000)
            plt.ylim(0,1000)
            plt.savefig('calibration_figures/corrected_uncalibrated_EvL_%s.png'%(i))
            i += 1
        plt.clf()

    

    def plot_uncorrected_calibrated_alpha_distributions(self):
        alphas, scale_factor = self.calibrate_uncorrected_alphas()
        tpcs = alphas.keys()
        color_codes = self.load_colors()
        colors = [color_codes[i] for i in range(0,len(tpcs))]
        i = 0
        for tpc in tpcs:
            plt.hist(alphas[tpc]['track_energy']/alphas[tpc]['length']*1e4, bins = 30, range = (200,800), linewidth = 2, edgecolor = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=1), color = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=0.2), label = tpc)
            i += 1
        plt.xlabel('dE/dx [keV/cm]')
        plt.legend()
        plt.title('Uncorrected Calibrated Horizontal Alpha Distribution')
        plt.savefig('calibration_figures/uncorrected_calibrated_alphas.png')
        plt.clf()
        
    def plot_corrected_calibrated_alpha_distributions(self):
        alphas, scale_factor = self.calibrate_corrected_alphas()
        tpcs = alphas.keys()
        color_codes = self.load_colors()
        colors = [color_codes[i] for i in range(0,len(tpcs))]
        i = 0
        for tpc in tpcs:
            plt.hist(alphas[tpc]['full_corrected_energy']/alphas[tpc]['length']*1e4, bins = 30, range = (200,800), linewidth = 2, edgecolor = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=1), color = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=0.2), label = tpc)
            i += 1
        plt.xlabel('dE/dx [keV/cm]')
        plt.legend()
        plt.title('Corrected Calibrated Horizontal Alpha Distribution')
        plt.savefig('calibration_figures/corrected_calibrated_alphas.png')
        plt.clf()

    def plot_corrected_calibrated_recoil_distributions(self):
        p = np.array([ 2.43367743e-06,  1.01652514e-02, -3.24232712e+00])
        recoils = self.calibrate_corrected_recoils()
        tpcs = recoils.keys()
        color_codes = self.load_colors()
        colors = [color_codes[i] for i in range(0,len(tpcs))]
        i = 0
        for tpc in tpcs:
            plt.plot(recoils[tpc]['length'], recoils[tpc]['full_corrected_energy'], 'o', markersize = 1, alpha = 0.5, label = tpc)
            plt.legend()
            plt.ylabel('Energy [keV]')
            plt.xlabel('length [um]')
            plt.xlim(0,20000)
            plt.ylim(0,1000)
            x = np.linspace(0,30000,10001)
            plt.savefig('calibration_figures/corrected_calibrated_EvL_%s.png'%(i))
            i += 1
        plt.plot(x, p[0]*x**2+p[1]*x+p[2])
        plt.savefig('calibration_figures/corrected_calibrated_EvL_all.png')
        plt.clf()
        i= 0
        recoils_uncorrected = self.calibrate_uncorrected_recoils()
        for tpc in tpcs: #plot uncorrected energy for reference
            plt.plot(recoils_uncorrected[tpc]['length'], recoils_uncorrected[tpc]['track_energy'], 'o', markersize = 1, alpha = 0.5, label = tpc)
            plt.legend()
            plt.ylabel('Energy [keV]')
            plt.xlabel('length [um]')
            plt.xlim(0,20000)
            plt.ylim(0,1000)
            plt.savefig('calibration_figures/uncorrected_calibrated_EvL_%s.png'%(i))
            i += 1
        plt.plot(x, p[0]*x**2+p[1]*x+p[2])
        plt.savefig('calibration_figures/uncorrected_calibrated_EvL_all.png')
        plt.clf()

    def post_calibration_cuts(self, dataframe):
        p = np.array([ 2.43367743e-06,  1.01652514e-02, -3.24232712e+00])
        dataframe = dataframe.loc[(dataframe['hitside_row_min']==0) & (dataframe['hitside_row_max']==0) & (dataframe['hitside_col_min']==0) & (dataframe['hitside_col_min']==0) & (dataframe['full_corrected_energy']>15) & (dataframe['full_corrected_energy']>(p[0]*dataframe['length']**2+p[1]*dataframe['length']+p[2]))]
        return dataframe

    def plot_corrected_calibrated_recoil_distributions_after_cuts(self):
        recoils = self.calibrate_corrected_recoils()
        tpcs = recoils.keys()
        for tpc in tpcs:
            recoils[tpc] = self.post_calibration_cuts(recoils[tpc])
        color_codes = self.load_colors()
        colors = [color_codes[i] for i in range(0,len(tpcs))]
        i = 0
        for tpc in tpcs:
            plt.plot(recoils[tpc]['length'], recoils[tpc]['full_corrected_energy'], 'o', markersize = 1, alpha = 0.5, label = tpc)
            plt.legend()
            plt.ylabel('Energy [keV]')
            plt.xlabel('length [um]')
            plt.xlim(0,10000)
            plt.ylim(0,500)
            plt.savefig('calibration_figures/corrected_calibrated_EvL_final_cuts_%s.png'%(i))
            i += 1
        plt.clf()
        i= 0
        for tpc in tpcs: #plot uncorrected energy for reference
            plt.plot(recoils[tpc]['length'], recoils[tpc]['full_corrected_energy'], 'o', markersize = 1, alpha = 0.5, label = tpc)
            plt.legend()
            plt.ylabel('Energy [keV]')
            plt.xlabel('length [um]')
            plt.xlim(0,10000)
            plt.ylim(0,500)
            plt.savefig('calibration_figures/uncorrected_calibrated_EvL_final_cuts_%s.png'%(i))
            i += 1
        plt.clf()

    def apply_recoil_cuts(self, corrected_energy = True, corrected_length = 0): #Cuts to train double Gaussian fit
        recoils =  self.calibrate_recoils(corrected_energy, corrected_length)
        for tpc in recoils.keys():
            if corrected_energy == False:
                #cut_min = np.array([ 5.51204819e-06, -6.90060241e-03,  1.78012048e+01])
                #cut_max = np.array([ 4.70883534e-06,  3.64728916e-02, -3.56124498e+01])
                cut_min = np.array([ 2.46987952e-06,  1.71265060e-02, -1.82530120e+01])
                cut_max = np.array([ 3.10240964e-06,  4.32198795e-02, -4.24397590e+01])
                ekey = 'track_energy'
                #recoils[tpc] = recoils[tpc].loc[recoils[tpc][ekey]>20]
                recoils[tpc] = recoils[tpc].loc[recoils[tpc][ekey]>8]
            else:
                #cut_min = np.array([ 4.26706827e-06,  1.33283133e-02, -9.98995984e+00])
                #cut_max = np.array([ 5.38152610e-06,  2.73975904e-02, -1.21285141e+01])
                if tpc == 'humu':
                    cut_min = np.array([1.62248996e-06, 3.18554217e-03, 5.89558233e+00])
                else:
                    cut_min = np.array([4.71887550e-06, 7.68072289e-03, 1.30522088e+00])
                cut_max = np.array([6.15461847e-06,  3.04006024e-02, -1.94678715e+01])
                ekey = 'full_corrected_energy'
                #recoils[tpc] = recoils[tpc].loc[recoils[tpc][ekey]>20]
                recoils[tpc] = recoils[tpc].loc[recoils[tpc][ekey]>10]
            index = recoils[tpc].loc[(recoils[tpc][ekey] > (cut_min[0]*recoils[tpc]['length']**2 + cut_min[1]*recoils[tpc]['length']+cut_min[2])) & (recoils[tpc][ekey]<(cut_max[0]*recoils[tpc]['length']**2  +cut_max[1]*recoils[tpc]['length']+cut_max[2]))].index.to_numpy()
            recoils[tpc] = recoils[tpc].loc[recoils[tpc][ekey] > (cut_min[0]*recoils[tpc]['length']**2 + cut_min[1]*recoils[tpc]['length']+cut_min[2])]
            recoils[tpc]['He_recoil'] = 0
            recoils[tpc]['He_recoil'][index] = 1
        return recoils

#legend_elements = [Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[0], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[0], alpha=1), label = 'iiwi'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[1], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[1], alpha=1), label = 'humu'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[2], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[2], alpha=1), label = 'nene'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[3], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[3], alpha=1), label = 'tako'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[4], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[4], alpha=1), label = 'palila'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[5], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[5], alpha=1), label = 'elepaio'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[0], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[0], alpha=1), hatch = '//\\\\', label = 'iiwi'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[1], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[1], alpha=1), hatch = '//\\\\', label = 'humu'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[2], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[2], alpha=1), hatch = '//\\\\', label = 'nene'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[3], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[3], alpha=1), hatch = '//\\\\', label = 'tako'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[4], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[4], alpha=1), hatch = '//\\\\', label = 'palila'), Patch(facecolor=matplotlib.colors.colorConverter.to_rgba(colors[5], alpha=0.2), edgecolor=matplotlib.colors.colorConverter.to_rgba(colors[5], alpha=1), hatch = '//\\\\', label = 'elepaio')]

'''
