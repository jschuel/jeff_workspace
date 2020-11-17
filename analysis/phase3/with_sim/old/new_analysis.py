import pandas as pd
import uproot as ur
import root_pandas as rp
import numpy as np
import matplotlib.pyplot as plt
import ROOT
from ROOT import TVector3
import array
from sklearn import linear_model
from sklearn.linear_model.base import LinearModel
from sklearn.base import RegressorMixin
from sklearn.utils import check_X_y
import numpy as np
import matplotlib
from matplotlib.lines import Line2D
from matplotlib import rc
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch

rc('text', usetex=False)

class analysis:

    def __init__(self, E_cut = 9, E_cut_err = 0, input_file= "~/data/phase3/spring_2020/05-09-20/combined_ntuples/05-09_whole_study_even_newerest.root", recoils_only = True, fei4_restrict = True): #enter negative value for E_Cut_err to get low systematic
        self.tpc_data = self.get_tpc_data(recoils_only = recoils_only, E_cut = E_cut, E_cut_err = E_cut_err)
        self.study_data = self.get_raw_study_data()
        self.MC_base = self.apply_chip_calibrations_to_MC(E_cut = E_cut)
        self.MC_data = self.create_VRCs()
        self.MC_rates = self.get_MC_rates(fei4_restrict = fei4_restrict, recoils_only = recoils_only)
        
    def get_tpc_data(self, input_dir = '~/data/phase3/spring_2020/05-09-20/tpc_root_files/', recoils_only = False, E_cut = 9, E_cut_err = 0):
        data = {}
        tpcs = ['iiwi', 'humu', 'nene', 'tako', 'palila', 'elepaio']
        for tpc in tpcs:
            data[tpc] = ur.open(input_dir + "%s_all_newester3.root"%(tpc))[ur.open(input_dir + "%s_all_newester3.root"%(tpc)).keys()[0]].pandas.df(flatten=False)
            data[tpc] = data[tpc].loc[data[tpc]['track_energy']>=(E_cut+E_cut_err)]
            if recoils_only == True:
                data[tpc] = data[tpc].loc[data[tpc]['is_recoil'] == 1]
        return data
    
    def get_raw_study_data(self, input_file= "~/data/phase3/spring_2020/05-09-20/combined_ntuples/05-09_whole_study_even_newerest.root"):
        study_data = ur.open(input_file)[ur.open(input_file).keys()[0]].pandas.df(flatten=False)
        #dfs = self.get_tpc_data()
        dfs = self.tpc_data
        tpcs = dfs.keys()
        for tpc in tpcs:
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

    def get_MC_data(self):
        tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all', 'RBB_Lumi', 'twoPhoton_Lumi']
        tree = 'tree_fe4_after_threshold'
        data = {}
        truth = {}
        for tpc in tpcs:
            data[tpc] = {}
            truth[tpc] = {}
            dir = '/home/jeef/data/phase3/spring_2020/05-09-20/geant4_simulation/all_events/%s/'%(tpc)
            #if bigChip == True and recoils_only == True:
            #    dir = '/home/jeef/data/phase3/spring_2020/05-09-20/geant4_simulation/big_chip/%s/'%(tpc)
            #elif bigChip == False and recoils_only == True:
            #    dir = '/home/jeef/data/phase3/spring_2020/05-09-20/geant4_simulation/%s/'%(tpc)
            #else:
            #    dir = '/home/jeef/data/phase3/spring_2020/05-09-20/geant4_simulation/all_events/%s/'%(tpc)
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
                    data[tpc][bg]['PDG'] = truth[tpc][bg]['tkPDG']
                    #if recoils_only == True:
                    #    data[tpc][bg]['truth_energy'] = truth[tpc][bg]['truthRecoilEnergy']
                    #    data[tpc][bg]['truth_mother_energy'] = truth[tpc][bg]['truthMotherEnergy']
                    #    data[tpc][bg][['truth_vertex_X', 'truth_vertex_Y', 'truth_vertex_Z', 'truth_px', 'truth_py', 'truth_pz', 'truth_mother_X', 'truth_mother_Y', 'truth_mother_Z', 'truth_mother_px', 'truth_mother_py', 'truth_mother_pz']] = truth[tpc][bg][['truthRecoilVtx_x_belle_frame', 'truthRecoilVtx_y_belle_frame', 'truthRecoilVtx_z_belle_frame', 'truthRecoilMom_x_belle_frame', 'truthRecoilMom_y_belle_frame', 'truthRecoilMom_z_belle_frame', 'truthMotherVtx_x_belle_frame', 'truthMotherVtx_y_belle_frame', 'truthMotherVtx_z_belle_frame', 'truthMotherMom_x_belle_frame', 'truthMotherMom_y_belle_frame', 'truthMotherMom_z_belle_frame']]
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

    def apply_chip_calibrations_to_MC(self, E_cut):
        MC = self.get_MC_data()
        tpcs = MC.keys()
        gain = {'iiwi': 1502, 'nene': 899, 'humu': 878, 'palila': 1033, 'tako': 807, 'elepaio': 797}
        W = 35.075
        tot_to_q = {}
        for tpc in tpcs:
            tot_to_q[tpc] = pd.DataFrame()
            tot_to_q[tpc]['tot_code'] = [i for i in range(0,14)]
            
        tot_to_q['iiwi']['conversion'] = np.array([1833.00, 2345.17, 3017.33, 6001.54, 8891.71,
                                  11497.43, 14335.32, 18081.33, 22526.06, 27236.90,
                                  32056.16, 36955.09, 41874.75, 46794.40])

        tot_to_q['nene']['conversion'] = np.array([2083.56, 2482.24, 4126.52, 5621.03, 7920.43,
                                                   11667.35, 15117.97, 19489.23, 23211.63, 27483.98,
                                                   32272.73, 37262.83, 42283.59, 47304.34])

        tot_to_q['humu']['conversion'] = np.array([2083.47, 2324.41, 3679.37, 5433.43, 6862.72,
                                                   10000.83, 13701.08, 17258.86, 21438.70, 25821.34,
                                                   30153.82, 34460.74, 39042.80, 43624.85])

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
            MC[tpc]['reco_energy'] = MC[tpc]['sumtot']*W/gain[tpc]*1e-3
            MC[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>=(E_cut)]

        return MC

    def determine_new_angles(self):
        def fitsvd_numpy(df, i): #faster than root fit, so this is standard
            x = df.iloc[i]['x']
            y = df.iloc[i]['y']
            z = df.iloc[i]['z']
            q = df.iloc[i]['q_from_tot']
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
        
        def update_ntuple(df):
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

        MC = self.MC_base
        tpcs = MC.keys()
        for tpc in tpcs:
            MC[tpc] = update_ntuple(MC[tpc])
        return MC

    def identify_MC_recoils(self):
        MC = self.determine_new_angles()
        x = {'iiwi': np.array([850, 2500, 15000]), 'nene': np.array([950, 2500, 16000]), 'humu': np.array([1950, 3000, 20000]),
             'palila': np.array([900, 2350, 15000]), 'tako': np.array([700, 2500, 15000]), 'elepaio': np.array([1050, 2500, 15000])}
        y = np.array([8,20,800])
        cut = {}
        tpcs = MC.keys()
        for tpc in tpcs:
            MC[tpc]['is_recoil'] = 0
            cut[tpc] = np.polyfit(x[tpc],y,2)
            index = MC[tpc].loc[(MC[tpc]['reco_energy']>=cut[tpc][0]*((MC[tpc]['fit_length']*10000)**2 + cut[tpc][1]*(MC[tpc]['fit_length']*10000) + cut[tpc][2]))].index.to_numpy()
            MC[tpc]['is_recoil'][index] = 1
        return MC

    def create_VRCs(self): #VRC=virtual readout chip. We create them when we have MC chip larger than fiducial area of chip
        
        def identify_edges(df):
            df['hit_edge'] = 0
            for i in range(0,len(df)):
                if ((0 in df['column'].iloc[i]) == True or (79 in df['column'].iloc[i]) == True 
                or (-78 in df['column'].iloc[i]) == True or (-159 in df['column'].iloc[i]) == True
                or (159 in df['column'].iloc[i]) == True or (0 in df['row'].iloc[i]) == True
                or (335 in df['row'].iloc[i]) == True or (-335 in df['row'].iloc[i]) == True
                or (670 in df['row'].iloc[i]) == True or (-670 in df['row'].iloc[i]) == True):
                    df['hit_edge'].iloc[i] = 1
            return df

        def Assign_VRC_ID(df):
            df['col_id'] = 0 #allowed values of 0-3 for 4x4 array of VRC
            df['row_id'] = 0 #allowed values of 0-3 for 4x4 array of VRC
            df['VRC_id'] = 0 #in practice we index from 1-9 given the locations of all hits
    
            index_col1 = df.loc[(df['column'].apply(lambda x: x.max()) > -79) & (df['column'].apply(lambda x: x.max()) <= 0)].index.to_numpy()
            index_col2 = df.loc[(df['column'].apply(lambda x: x.max()) > 0) & (df['column'].apply(lambda x: x.max()) <= 79)].index.to_numpy()
            index_col3 = df.loc[(df['column'].apply(lambda x: x.max()) > 79) & (df['column'].apply(lambda x: x.max()) <= 158)].index.to_numpy()
            df['col_id'][index_col1] = 1
            df['col_id'][index_col2] = 2
            df['col_id'][index_col3] = 3
    
            index_row1 = df.loc[(df['row'].apply(lambda x: x.max()) > -335) & (df['row'].apply(lambda x: x.max()) <= 0)].index.to_numpy()
            index_row2 = df.loc[(df['row'].apply(lambda x: x.max()) > 0) & (df['row'].apply(lambda x: x.max()) <= 335)].index.to_numpy()
            index_row3 = df.loc[(df['row'].apply(lambda x: x.max()) > 335) & (df['row'].apply(lambda x: x.max()) <= 670)].index.to_numpy()
            df['row_id'][index_row1] = 1
            df['row_id'][index_row2] = 2
            df['row_id'][index_row3] = 3
            for i in range(1,4):
                index_top = df.loc[(df['col_id']==i) & (df['row_id']==3)].index.to_numpy()
                df['VRC_id'][index_top] = int(i)
                index_mid = df.loc[(df['col_id']==i) & (df['row_id']==2)].index.to_numpy()
                df['VRC_id'][index_mid] = i+3
                index_bot = df.loc[(df['col_id']==i) & (df['row_id']==1)].index.to_numpy()
                df['VRC_id'][index_bot] = i+6

            return df

        MC = self.identify_MC_recoils()
        tpcs = MC.keys()
        MC_red = {}
        for tpc in tpcs:
            MC[tpc] = identify_edges(MC[tpc])
            MC[tpc]['ones'] = 1
            MC_red[tpc] = MC[tpc].loc[MC[tpc]['hit_edge'] == 0] #edge cut
            MC_red[tpc].index = [i for i in range(0,len(MC_red[tpc]))]
            MC_red[tpc] = Assign_VRC_ID(MC_red[tpc])

        return MC_red
    
    def get_MC_rates(self, fei4_restrict = True, recoils_only = True, I_HER = 1000, I_LER = 1200, sy_LER=37, sy_HER=36, nb_LER=1576, nb_HER=1576, lumi=25): #Scale to luminosity of interest. Units: 1e34cm-2s-1
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
            #MC[tpc] = MC[tpc].loc[MC[tpc]['reco_energy']>E_cut]
            if recoils_only == True and fei4_restrict == True:
                MC_new[tpc] = MC[tpc].loc[MC[tpc]['is_recoil'] == 1]
                MC_new[tpc] = MC_new[tpc].loc[MC_new[tpc]['VRC_id']==5]
            elif recoils_only == False and fei4_restrict == True:
                MC_new[tpc] = MC[tpc].loc[MC[tpc]['VRC_id']==5]
            elif recoils_only == True and fei4_restrict == False:
                MC_new[tpc] = MC[tpc].loc[MC[tpc]['is_recoil'] == 1]
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
                    t = 1e-3*lumi_frac                    
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
        df['LER_T'] = df['Touschek_LER_all']/(1200/(37*1576))*(I_LER/(sy_LER*nb_LER))
        df['HER_T'] = df['Touschek_HER_all']/(1000/(36*1576))*(I_HER/(sy_HER*nb_HER))
        df['Lumi'] = df['RBB_Lumi'] + df['twoPhoton_Lumi']
        df = df[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
        df_err.index = tpcs
        df_err['LER_bg_base'] = (df_err['Brems_LER_base'] + df_err['Coulomb_LER_base'])/1200*I_LER #scale by input I/I_ref
        df_err['HER_bg_base'] = (df_err['Brems_HER_base'] + df_err['Coulomb_HER_base'])/1000*I_HER
        df_err['LER_bg_dynamic'] = (df_err['Brems_LER_dynamic'] + df_err['Coulomb_LER_dynamic'])/1200**2*I_LER**2
        df_err['HER_bg_dynamic'] = (df_err['Brems_HER_dynamic'] + df_err['Coulomb_HER_dynamic'])/1000**2*I_HER**2
        df_err['LER_T'] = df_err['Touschek_LER_all']/(1200/(37*1576))*(I_LER/(sy_LER*nb_LER))
        df_err['HER_T'] = df_err['Touschek_HER_all']/(1000/(36*1576))*(I_HER/(sy_HER*nb_HER))
        df_err['Lumi'] = np.sqrt(df_err['RBB_Lumi']**2 + df_err['twoPhoton_Lumi']**2)
        df_err = df_err[['LER_bg_base', 'LER_bg_dynamic', 'LER_T', 'HER_bg_base', 'HER_bg_dynamic', 'HER_T', 'Lumi']]
        df_err.columns = ['LER_bg_base_err', 'LER_bg_dynamic_err', 'LER_T_err', 'HER_bg_base_err', 'HER_bg_dynamic_err', 'HER_T_err', 'Lumi_err']
        df_combined = pd.concat([df,df_err], axis=1)
        return df_combined

    def select_study(self, study_type, study_period): #LER, HER, Lumi, Cont_inj, Decay
        #raw_data = self.get_raw_study_data()
        raw_data = self.study_data
        study_data = raw_data.loc[(raw_data['%s_study_flag'%(study_type)]==1) & (raw_data['%s_flag'%(study_period)] == 1)]
        return study_data

    def get_tpc_data_during_study_period(self, study_type, study_period):
        study_data = self.select_study(study_type, study_period)
        tpc_data = self.tpc_data
        tpcs = tpc_data.keys()
        tpc_study_data = {}
        for tpc in tpcs:
            tpc_study_data[tpc] = tpc_data[tpc].loc[(tpc_data[tpc]['timestamp_start'].astype('int')).isin(study_data['ts'].astype('int'))==True]
        return tpc_study_data

    def partition_data_into_subsets(self, study_type, study_period, bins=12):
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

    def compute_means_and_errs(self, study_type, study_period, bins=12):
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

    def get_fit_parameters(self, study_type, study_period, bins=12): #Gives parameters B0, B1, and T defined by Rate/I = B0 + B1*I + T*I/(sy*Nb)
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
            #f2.SetParLimits(0,2e-12,2e-4)
            #f2.SetParLimits(1,2e-9,2e-2)
            #f2.SetParLimits(2,1e-12,1)
            #f2.SetParLimits(0,1e-12,1e-6)
            #f2.SetParLimits(1,1e-9,1e-2)
            #f2.SetParLimits(2,1e-13,2e-1)
            #f2.SetParLimits(0,1e-12,1e-5)
            #f2.SetParLimits(1,1e-8,1e-2)
            #f2.SetParLimits(2,1e-8,1)
            f2.SetParLimits(0,0,1e-5)
            f2.SetParLimits(1,0,1e-2)
            f2.SetParLimits(2,0,1)
            gr = ROOT.TGraph2DErrors(len(x1), x1, x2, y_root, x1err, x2err, y_root_err)
            gr.Fit(f2, 'SREM')
            fit[tpc+'_B0'] = f2.GetParameter(0)
            fit[tpc+'_B1'] = f2.GetParameter(1)
            fit[tpc+'_T'] = f2.GetParameter(2)
            fit[tpc+'_B0_err'] = f2.GetParError(0)
            fit[tpc+'_B1_err'] = f2.GetParError(1)
            fit[tpc+'_T_err'] = f2.GetParError(2)
        return fit

    def plot_fit(self, tpc, study_type,  bins=12, ymax = 2, legend = False):
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
        #plt.ylim(-0.01,ymax)
        plt.ylabel('Rate[Hz]')
        if legend == True:
            #if study_type == "LER":
            #    custom_lines = [Line2D([0], [0], color='red', lw=4, label = '%s Current [mA]'%(study_type)), Line2D([0], [0], color='magenta', lw=4, label = 'TPC predicted rate'), Line2D([0], [0], color='limegreen', lw=4, label = 'TPC predicted Touschek rate'), Line2D([0], [0], color='cyan', lw=4, label = 'TPC predicted beam-gas rate'), Line2D([0], [0], marker='o', color='w', markersize=6, lw=0, markerfacecolor = 'black', label = 'Data')]
            #else:
            #    custom_lines = [Line2D([0], [0], color='blue', lw=4, label = '%s Current [mA]'%(study_type)), Line2D([0], [0], color='magenta', lw=4, label = 'TPC predicted rate'), Line2D([0], [0], color='limegreen', lw=4, label = 'TPC predicted Touschek rate'), Line2D([0], [0], color='cyan', lw=4, label = 'TPC predicted beam-gas rate'), Line2D([0], [0], marker='o', color='w', markersize=6, lw=0, markerfacecolor = 'black', label = 'Data')]
            #plt.legend(handles = custom_lines, loc = 'upper right')
            plt.legend([p1,p2,p3,p4,p5],['Data', 'Total Fit', 'Beam Gas Base', 'Beam Gas Dyn.', 'Touschek'])
    def make_fwd_plots(self, study_type, bins=12):
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
        plt.subplot(3,1,1)
        self.plot_fit("iiwi", study_type, bins, ymax=1, legend = True)
        plt.title("Iiwi (z = +6.6 m)")
        plt.xticks([])
        plt.subplot(3,1,2)
        self.plot_fit("nene", study_type, bins, ymax)
        plt.title("Nene (z = +14 m)")
        plt.xticks([])
        plt.subplot(3,1,3)
        self.plot_fit("humu", study_type, bins, ymax)
        plt.title("Humu (z = +16 m)")
        plt.xticks([])
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0.2)
        plt.subplots_adjust(wspace = 0.5)
        plt.show()

    def make_bwd_plots(self, study_type, bins=12):
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

    def measure_and_fit_lumi_bgs(self, tpc, study_period ,bins=12):
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

    def plot_all_luminosity(self, study_period = "Cont_inj", bins=12):
        plt.figure(figsize = (15,15))
        plt.rc('legend', fontsize=12)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=16)
        plt.rc('axes', titlesize=16)
        plt.subplot(3,2,1)
        self.measure_and_fit_lumi_bgs("palila", study_period, bins)
        plt.title("z = -5.6m")
        plt.grid()
        plt.subplot(3,2,2)
        self.measure_and_fit_lumi_bgs("iiwi", study_period, bins)
        plt.title("z = +6.6m")
        plt.grid()
        plt.subplot(3,2,3)
        self.measure_and_fit_lumi_bgs("tako", study_period, bins)
        plt.title("z = -8.0m")
        plt.grid()
        plt.subplot(3,2,4)
        self.measure_and_fit_lumi_bgs("nene", study_period, bins)
        plt.title("z = +14m")
        plt.grid()
        plt.subplot(3,2,5)
        self.measure_and_fit_lumi_bgs("elepaio", study_period, bins)
        plt.title("z = -14m")
        plt.grid()
        plt.subplot(3,2,6)
        self.measure_and_fit_lumi_bgs("humu", study_period, bins)
        plt.title("z = +16m")
        plt.subplots_adjust(hspace=0.5)
        plt.subplots_adjust(wspace=0.5)
        plt.grid()
        #plt.savefig('lumi_fits.png', bbox_inches='tight')
        plt.show()

    def plot_bwd_luminosity(self, bins=12):
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

    def plot_fwd_luminosity(self, bins=12):
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

    def plot_bg_summary(self, study_period, bins=12, MC = False, I_HER = 1000, I_LER = 1200, sy_LER=37, sy_HER=36, nb_LER=1576, nb_HER=1576, L=25):
        #tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        if MC == True:
            #df = self.get_MC_rates(E_cut = 8, I_HER = I_HER, I_LER = I_LER, sy_LER=sy_LER, sy_HER=sy_HER, nb_LER=nb_LER, nb_HER=nb_HER, lumi=L)
            df = self.MC_rates
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
            for tpc in tpcs:
                #I_HER = 1000
                #I_LER = 1200
                #sy_LER = 37
                #sy_HER = 36
                #nb_LER = 1576
                #nb_HER = 1576
                #L = 25 #luminoisty in 1e34 cm-2s-1
                lumi_fits = self.measure_and_fit_lumi_bgs(tpc, study_period ,bins=20)
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

    def compute_data_MC_ratios(self, study_period, bins=12, E_cut = 8, I_HER = 1000, I_LER = 1200, sy_LER=37, sy_HER=36, nb_LER=1576, nb_HER=1576, L=25):

        #MC = self.get_MC_rates(E_cut = E_cut, I_HER = I_HER, I_LER = I_LER, sy_LER=sy_LER, sy_HER=sy_HER, nb_LER=nb_LER, nb_HER=nb_HER, lumi=L)
        MC = self.MC_rates
        tpcs = ['elepaio', 'tako', 'palila', 'iiwi', 'nene', 'humu']
        LER_fit_params = self.get_fit_parameters("LER", study_period, bins)
        HER_fit_params = self.get_fit_parameters("HER", study_period, bins)
        fit_dict = {}
        df = pd.DataFrame() #order is LER_bg, LER_T, HER_bg, HER_T, Lumi 
        for tpc in tpcs:
            #I_HER = 1000
            #I_LER = 1200
            #sy_LER = 37
            #sy_HER = 36
            #nb_LER = 1576
            #nb_HER = 1576
            #L = 25 #luminoisty in 1e34 cm-2s-1
            lumi_fits = self.measure_and_fit_lumi_bgs(tpc, study_period ,bins=20)
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

a = analysis()
