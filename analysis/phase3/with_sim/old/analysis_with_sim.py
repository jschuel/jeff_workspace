import uproot as ur
import ROOT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

class analysis:

    def __init__(self):
        #self.data = self.get_MC_data()
        #self.rates = self.get_simulated_rates()
        pass
    def get_simulated_rates(self):
        tpcs = ['iiwi', 'palila', 'tako', 'elepaio', 'nene', 'humu']
        bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all', 'RBB_Lumi', 'twoPhoton_Lumi']
        tree = 'tree_fe4_after_threshold'
        rates = {}
        df = pd.DataFrame()
        for tpc in tpcs:
            rates[tpc] = {}
            dir = '/home/jeef/data/phase3/spring_2020/05-09-20/geant4_simulation/%s/'%(tpc)
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
                    t = 6e-5
                elif (bg == 'twoPhoton_Lumi'):
                    t = 0.005
                try:
                    rates[tpc][bg] = len(ur.open(dir+bg+'_%s.root'%(tpc))[tree].pandas.df(flatten=False))/(t*100.) #100 compensates for modified cross section
                except FileNotFoundError:
                    rates[tpc][bg] = 0
            df = df.append(pd.DataFrame.from_dict(rates[tpc], 'index').T)
        df.index = tpcs
        return df

    def get_MC_data(self):
        #tpcs = ['iiwi', 'palila', 'tako', 'elepaio', 'nene', 'humu']
        tpcs = ['palila']
        bgtype = ['RBB_Lumi']
        #bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all', 'RBB_Lumi', 'twoPhoton_Lumi']
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
            dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('Coulomb', 'BeamGas', regex=False)
            dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('Brems', 'BeamGas', regex=False)
            dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('_all', '', regex=False)
            dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('twoPhoton_', '', regex=False)
            dfs[tpc]['bgType'] = dfs[tpc]['bgType'].str.replace('RBB_', '', regex=False)
            
        return dfs

    def visualize_MC_with_geometry(self):
        plt.rc('legend', fontsize=12)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=18)
        plt.rc('axes', titlesize=18)

        cm = matplotlib.cm.plasma
        norm = matplotlib.colors.LogNorm()
        #norm = matplotlib.colors.Normalize(vmin = 0)
        sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)

        data = self.get_MC_data()[0]
        
        for key in data.keys():
            data[key] = data[key].loc[(data[key]['truth_energy']>0) & (data[key]['truth_energy'].duplicated() == False)]
            data[key].index = [i for i in range(0,len(data[key]))]

        ###Make Consistent Color Scale###
        df = pd.DataFrame()
        for key in data.keys():
            df = df.append(data[key])
        df.index = [i for i in range(0,len(df))]
        color = cm(norm(df['truth_mother_energy']))
        #color = cm(norm(df['truth_energy']/df['truth_mother_energy']))
        ###
        
        img = plt.imread("/home/jeef/Pictures/farbeamline_nocolor.png")
        fig, ax = plt.subplots(2,1,figsize = (20,16))
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        ax[1].set_xticks([])
        ax[1].set_yticks([])
        #ax[0].imshow(np.flipud(img), origin = 'lower', extent = [-3300,3170, -267,233], aspect = 'auto')
        #ax[1].imshow(np.flipud(img), origin = 'lower', extent = [-3300,3170, -267,233], aspect = 'auto')
        ax[0].imshow(np.flipud(img), origin = 'lower', extent = [-3333,3142, -438,414], aspect = 'auto')
        ax[1].imshow(np.flipud(img), origin = 'lower', extent = [-3333,3142, -438,414], aspect = 'auto')
        j = 0
        for key in data.keys():
            ax[0].plot(data[key]['truth_vertex_Z'].mean(), data[key]['truth_vertex_X'].mean(), 's', color = 'green', markersize = 20, alpha = 0.5)
            ax[1].plot(data[key]['truth_vertex_Z'].mean(), data[key]['truth_vertex_X'].mean(), 's', color = 'green', markersize = 20, alpha = 0.5)
            for i in range(0,len(data[key])):
                ax[0].scatter(data[key]['truth_mother_Z'].iloc[i], data[key]['truth_mother_X'].iloc[i], c = color[j])
                zs = [data[key]['truth_mother_Z'].iloc[i], data[key]['truth_vertex_Z'].iloc[i]]
                xs = [data[key]['truth_mother_X'].iloc[i], data[key]['truth_vertex_X'].iloc[i]]
                ax[1].plot(zs, xs, color = color[j])
                j+=1
        
        plt.colorbar(sm, ax=ax[0]).set_label('Neutron Energy [keV]', rotation = 270, labelpad = 20)
        plt.colorbar(sm, ax=ax[1]).set_label('Neutron Energy [keV]', rotation = 270, labelpad = 20)
        #plt.colorbar(sm, ax=ax[0]).set_label(r'$E_{recoil}/E_{neutron}$', rotation = 270, labelpad = 20)
        #plt.colorbar(sm, ax=ax[1]).set_label(r'$E_{recoil}/E_{neutron}$', rotation = 270, labelpad = 20)
        
        #plt.savefig("neutron_production_points.png")
        #plt.savefig("neutron_production_points_ratio.png")
        plt.show()

    def visualize_3D(self):
        cm = matplotlib.cm.plasma
        norm = matplotlib.colors.LogNorm(vmin = 100)
        #norm = matplotlib.colors.Normalize(vmin = 0)
        sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)

        data = self.get_MC_data()[0]
        
        for key in data.keys():
            data[key] = data[key].loc[(data[key]['truth_energy']>0) & (data[key]['truth_energy'].duplicated() == False)]
            data[key].index = [i for i in range(0,len(data[key]))]

        ###Make Consistent Color Scale###
        df = pd.DataFrame()
        for key in data.keys():
            df = df.append(data[key])
        df.index = [i for i in range(0,len(df))]
        color = cm(norm(df['truth_mother_energy']))
        #color = cm(norm(df['truth_energy']/df['truth_mother_energy']))
        ###

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        j=0
        for key in data.keys():
            ax.plot(data[key]['truth_vertex_Z'].mean(), data[key]['truth_vertex_X'].mean(), data[key]['truth_vertex_Y'].mean(),'s', color = 'green', markersize = 20)
            for i in range(0,len(data[key])):
                zs = [data[key]['truth_mother_Z'].iloc[i], data[key]['truth_vertex_Z'].iloc[i]]
                xs = [data[key]['truth_mother_X'].iloc[i], data[key]['truth_vertex_X'].iloc[i]]
                ys = [data[key]['truth_mother_Y'].iloc[i], data[key]['truth_vertex_Y'].iloc[i]]
                ax.plot(zs, xs, ys, color = color[j])
                j+=1
        #ax.scatter(df['truth_mother_Z'], df['truth_mother_X'], df['truth_mother_Y'], c = color)
        ax.set_xlabel('z')
        ax.set_ylabel('x')
        ax.set_zlabel('y')
        plt.colorbar(sm).set_label('Neutron Energy [keV]')
        plt.show()
        
a = analysis()
#data, truth = a.get_MC_data()
a.visualize_MC_with_geometry()
#a.visualize_3D()
df = a.get_simulated_rates()
