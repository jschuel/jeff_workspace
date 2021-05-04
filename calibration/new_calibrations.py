
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
            
class tpc_calibration:

    def __init__(self):
        self.updated_recoils, self.updated_alphas = self.determine_new_angles_and_lengths()
        self.calibrated_alphas, self.scale_factor, self.scale_factor_err = self.calibrate_alphas()
        self.calibrated_recoils = self.calibrate_recoils()
        self.plot_alpha_distributions(calibrated = True, gain = True)
        #self.write_to_root_file()
        #pass
        
    def get_tpc_list(self, tpc_list = ['iiwi', 'humu', 'nene', 'tako', 'palila', 'elepaio']):
        return tpc_list

    def base_alphas(self):
        tpcs = self.get_tpc_list()
        df = {}
        for tpc in tpcs:
            #df[tpc] = ur.open('~/data/phase3/spring_2020/maintenance_day_test/%s_alphas_all.root'%(tpc))[ur.open('~/data/phase3/spring_2020/maintenance_day_test/%s_alphas_all.root'%(tpc)).keys()[0]].pandas.df(flatten=False)
            df[tpc] = rp.read_root('~/data/phase3/spring_2020/maintenance_day_test/%s_alphas_all.root'%(tpc))
            df[tpc]['track_energy'] = df[tpc]['track_charge']*34.4525/2000*1e-3 # "uncalibrates" energy to start fresh
            df[tpc]['pixel_energy'] = df[tpc]['pixel_charge']*34.4525/2000*1e-3 # "uncalibrates" energy to start fresh
        return df #returns dictionary of horizontal alphas

    def get_study_data(self):
        tpcs = self.get_tpc_list()
        df = {}
        for tpc in tpcs:
            #df[tpc] = ur.open('~/data/phase3/spring_2020/05-09-20/tpc_root_files/%s_all_new.root'%(tpc))[ur.open('~/data/phase3/spring_2020/05-09-20/tpc_root_files/%s_all_new.root'%(tpc)).keys()[0]].pandas.df(flatten=False)
            df[tpc] = rp.read_root('~/data/phase3/spring_2020/05-09-20/tpc_root_files/%s_all_new.root'%(tpc))
            df[tpc]['track_energy'] = df[tpc]['track_charge']*34.4525/2000*1e-3 # "uncalibrates" energy to start fresh
            df[tpc]['pixel_energy'] = df[tpc]['pixel_charge']*34.4525/2000*1e-3 # "uncalibrates" energy to start fresh
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

    def determine_new_angles_and_lengths(self):
        recoils = self.load_recoils()
        alphas = self.base_alphas()
        tpcs = recoils.keys()
        v_drift = {'iiwi': 216.25, 'nene': 216.25, 'humu': 216.25, 'palila': 216.25, 'tako': 216.25, 'elepaio': 151.675}
        v_drift_alphas = {'iiwi': 216.25, 'nene': 216.25, 'humu': 216.25, 'palila': 216.25, 'tako': 259.25, 'elepaio': 173.4}
        ### Scale z by correct drift speed ###
        zs_alphas = {}
        for tpc in tpcs:
            zs_alphas[tpc] = []
            recoils[tpc]['z'] = recoils[tpc]['BCID']*v_drift[tpc]
            alphas[tpc]['z'] = alphas[tpc]['BCID']*v_drift_alphas[tpc]
            if tpc == 'elepaio':
                for row in alphas[tpc].itertuples():
                    if (row.timestamp_start > 1587421058) and (row.timestamp_start < 1587853058):
                        zs_alphas[tpc].append(row.BCID*164.925)
                    elif (row.timestamp_start > 1587853058):
                        zs_alphas[tpc].append(row.BCID*151.675)
                    else:
                        zs_alphas[tpc].append(row.BCID*v_drift_alphas[tpc])
                alphas[tpc]['z'] = zs_alphas[tpc]
            if tpc == 'tako':
                for row in alphas[tpc].itertuples():
                    if (row.timestamp_start > 1587853058):
                        zs_alphas[tpc].append(row.BCID*216.25)
                    else:
                        zs_alphas[tpc].append(row.BCID*v_drift_alphas[tpc])
                alphas[tpc]['z'] = zs_alphas[tpc]
                
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
            print("Performing SVD fits for %s RECOILS"%(tpc))
            recoils[tpc] = update_dataframe(recoils[tpc])
            print("DONE")
            print("Performing SVD fits for %s ALPHAS"%(tpc))
            alphas[tpc] = update_dataframe(alphas[tpc])
            print("DONE")
        return recoils, alphas

    def horizontal_alphas(self):
        alphas = self.updated_alphas
        df = {}
        tpcs = alphas.keys()
        for tpc in tpcs:
            df[tpc] = alphas[tpc].loc[(np.abs(alphas[tpc]['new_phi']*180/np.pi) < 5) | (np.abs((alphas[tpc]['new_phi'])*180/np.pi-180) < 5) | (np.abs(-180 - (alphas[tpc]['new_phi']*180/np.pi)) < 5)]
            df[tpc] = df[tpc].loc[(df[tpc]['new_theta']*180/np.pi>85) & (df[tpc]['new_theta']*180/np.pi<95)]
            #df[tpc]['BCID_range'] = [df[tpc]['BCID'][i].max()-df[tpc]['BCID'][i].min() for i in range (0,len(df[tpc]))]
            df[tpc] = df[tpc].loc[(df[tpc]['track_energy'] > 400)]# & (df[tpc]['BCID_range'] < 90)] #selection for alphas
            df[tpc] = df[tpc].loc[(np.abs(df[tpc]['phi']) < 5) | (np.abs(df[tpc]['phi']-180) < 5) | (np.abs(-180 - df[tpc]['phi']) < 5)]
            df[tpc] = df[tpc].loc[(df[tpc]['theta']>85) & (df[tpc]['theta']<95)]
            df[tpc]['track_energy'] = df[tpc]['track_charge']*34.4525/2000*1e-3 # "uncalibrates" energy to start fresh
            df[tpc]['pixel_energy'] = df[tpc]['pixel_charge']*34.4525/2000*1e-3 # "uncalibrates" energy to start fresh
        return df #returns dictionary of horizontal alphas


    def calibrate_alphas(self): #calibrates alphas to a de/dx reference value of 500 keV/cm
        alphas = self.horizontal_alphas()
        ekey = 'track_energy'
        tpcs = alphas.keys()
        scale_factor = {}
        scale_factor_err = {}
        for tpc in tpcs:
            #alphas[tpc]['dedx'] = alphas[tpc][ekey]/alphas[tpc][lkey]
            scale_factor[tpc] = 1430/((alphas[tpc][ekey]*np.sin(alphas[tpc]['new_theta'])).mean())
            scale_factor_err[tpc] = scale_factor[tpc]*(alphas[tpc][ekey]*np.sin(alphas[tpc]['new_theta'])).std()/((alphas[tpc][ekey]*np.sin(alphas[tpc]['new_theta'])).mean())
            #print((alphas[tpc][ekey]*np.sin(alphas[tpc]['theta']*np.pi/180)).sem())
            alphas[tpc][ekey] = scale_factor[tpc]*alphas[tpc][ekey]
            alphas[tpc][ekey+'_err'] = scale_factor[tpc]*alphas[tpc][ekey]*(alphas[tpc][ekey]*np.sin(alphas[tpc]['new_theta'])).std()/((alphas[tpc][ekey]*np.sin(alphas[tpc]['new_theta'])).mean())
            #print(tpc, scale_factor_err[tpc]/scale_factor[tpc])
        return alphas, scale_factor, scale_factor_err

    def calibrate_recoils(self): #calibrates recoils to a de/dx reference value of 500 keV/cm
        recoils = self.updated_recoils
        ekey = 'track_energy'
        scale_factor, scale_factor_err = self.scale_factor, self.scale_factor_err
        tpcs = recoils.keys()
        for tpc in tpcs:
            recoils[tpc]['track_energy'] = scale_factor[tpc]*recoils[tpc]['track_energy']
            recoils[tpc]['track_energy_err'] = scale_factor_err[tpc]*recoils[tpc]['track_energy']
            #recoils[tpc]['track_energy_err'] = scale_factor_err[tpc]*alphas[tpc]['track_energy'].mean()
            #recoils[tpc]['track_energy_err'] = scale_factor_err[tpc]*recoils[tpc]['track_energy']/scale_factor[tpc]
            #print(tpc, recoils[tpc]['track_energy_err']/recoils[tpc]['track_energy'])
        return recoils
            
    def load_colors(self): #makes colorblind friendly palette
        color_codes = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00'] #friendly colors for color blind
        #color_codes = ['#FFD700', '#4B0082', '#FF00FF', '#00FFFF', '#FF0000', '#808080']
        return color_codes

    def plot_alpha_distributions(self, calibrated = True, gain = False):
        if calibrated == True:
            alphas, scale_factor, scale_factor_err = self.calibrate_alphas()
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
                alphas[tpc][ekey] = alphas[tpc]['track_charge']*34.4525/alphas[tpc][ekey]*1e-3
                print(tpc, alphas[tpc][ekey].mean())
                #print(tpc, alphas[tpc][ekey].unique(), alphas[tpc][ekey+'_err'].unique())
                print(tpc, alphas[tpc][ekey].mean(), scale_factor_err[tpc]*alphas[tpc][ekey].mean())
                plt.bar(locations[i], alphas[tpc][ekey].mean(), yerr = scale_factor_err[tpc]*alphas[tpc][ekey].mean(), linewidth = 2, edgecolor = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=1), color = matplotlib.colors.colorConverter.to_rgba(colors[i], alpha=0.2))
                plt.xlabel('TPC Location')
                plt.ylabel(r'$G_{eff}$')
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
        
    def apply_recoil_cuts(self): #Cuts without xray veto
        recoils = self.calibrated_recoils
        #x = {'iiwi': np.array([850, 2500, 15000]), 'nene': np.array([950, 2500, 16000]), 'humu': np.array([1950, 3000, 20000]),
        #     'palila': np.array([900, 2350, 15000]), 'tako': np.array([1400, 2500, 15000]), 'elepaio': np.array([1050, 2500, 15000])}
        #y = np.array([6,20,800])
        #cut = {}
        xs = np.array([1200, 4000,10000])
        ys = np.array([7,50,280])
        cut = np.polyfit(xs,ys,2)
        tpcs = recoils.keys()
        for tpc in tpcs:
            #cut[tpc] = np.polyfit(x[tpc],y,2)
            recoils[tpc]['is_recoil'] = 0
            index = recoils[tpc].loc[(recoils[tpc]['track_energy']>=(cut[0]*recoils[tpc]['new_length']**2 + cut[1]*recoils[tpc]['new_length'] + cut[2]))].index.to_numpy()
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
            output = ROOT.TFile(outdir + '%s_all_newester7.root'%(tpc), 'recreate')
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
