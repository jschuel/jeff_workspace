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
        tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all']
        tree = 'tree_fe4_after_threshold'
        rates = {}
        df = pd.DataFrame()
        for tpc in tpcs:
            rates[tpc] = {}
            dir = '/home/jeef/data/phase3/spring_2020/05-09-20/geant4_simulation/truth_neutrons_only/%s/'%(tpc)
            for bg in bgtype:
                if (bg == 'Brems_HER_base') or (bg == 'Brems_HER_dynamic'):
                    t = 8.0
                elif (bg == 'Coulomb_HER_base') or (bg == 'Coulomb_HER_dynamic'):
                    t = 0.8
                elif (bg == 'Brems_LER_base') or (bg == 'Brems_LER_dynamic'):
                    t = 4.0
                elif (bg == 'Coulomb_LER_base') or (bg == 'Coulomb_LER_dynamic'):
                    t = 0.4
                elif (bg == 'Touschek_HER_all'):
                    t = 0.8
                elif (bg == 'Touschek_LER_all'):
                    t = 0.1
                try:
                    rates[tpc][bg] = len(ur.open(dir+bg+'_%s.root'%(tpc))[tree].pandas.df(flatten=False))/(t*100.) #100 compensates for modified cross section
                except FileNotFoundError:
                    rates[tpc][bg] = 0
            df = df.append(pd.DataFrame.from_dict(rates[tpc], 'index').T)
        df.index = tpcs
        return df

    def get_MC_data(self, bgType):
        tpcs = ['iiwi', 'palila', 'tako', 'elepaio']
        #bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all']
        #bgtype = ['twoPhoton_Lumi']
        truth = {}
        for tpc in tpcs:
            dir = '/home/jeef/data/phase3/spring_2020/05-09-20/geant4_simulation/truth_neutrons_only/%s/'%(tpc)
            try:
                try:
                    truth[tpc+'_'+bgtype] = ur.open(dir+bgtype+'_'+tpc+'_all.root')['recoils'].pandas.df(flatten=False)
                except KeyError:
                    truth[tpc+'_'+bgtype] = ur.open(dir+bgtype+'_'+tpc+'_all.root')['my_ttree'].pandas.df(flatten=False)
                try:
                    truth[tpc+'_'+bgtype]['chipx'] = truth[tpc+'_'+bgtype]['chipx'].apply(lambda x: x[0])
                    truth[tpc+'_'+bgtype]['chipy'] = truth[tpc+'_'+bgtype]['chipy'].apply(lambda x: x[0])
                    truth[tpc+'_'+bgtype]['chipz'] = truth[tpc+'_'+bgtype]['chipz'].apply(lambda x: x[0])
                except TypeError:
                    pass
            except FileNotFoundError:
                truth[tpc+'_'+bgtype] = pd.DataFrame()
        return truth

    def visualize_MC_with_geometry(self, bgtype):
        plt.rc('legend', fontsize=12)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=18)
        plt.rc('axes', titlesize=18)

        cm = matplotlib.cm.viridis
        norm = matplotlib.colors.LogNorm(vmin = 1e1, vmax = 1e5)
        #norm = matplotlib.colors.Normalize(vmin = 0)
        sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)

        data = self.get_MC_data(bgtype)
        print(data)
        
        ###Make Consistent Color Scale###
        df = pd.DataFrame()
        for key in data.keys():
            df = df.append(data[key])
        df.index = [i for i in range(0,len(df))]
        color = cm(norm(df['truthNeutronEnergy']))
        ###
        
        img = plt.imread("/home/jeef/Pictures/farbeamline.png")
        fig, ax = plt.subplots(2,1,figsize = (12,10))
        fig.suptitle('%s neutron SimHits and vertices'%(bgtype), fontsize=16)
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        ax[1].set_xticks([])
        ax[1].set_yticks([])
        ax[0].imshow(np.flipud(img), origin = 'lower', extent = [-3300,3170, -267,233], aspect = 'auto')
        ax[1].imshow(np.flipud(img), origin = 'lower', extent = [-3300,3170, -267,233], aspect = 'auto')
        j = 0
        for key in data.keys():
            for i in range(0,len(data[key])):
                ax[0].scatter(data[key]['truthNeutronVtx_z_belle_frame'].iloc[i], data[key]['truthNeutronVtx_x_belle_frame'].iloc[i], c = color[j], s=0.8)
                zs = [data[key]['truthNeutronVtx_z_belle_frame'].iloc[i], data[key]['chipz'].iloc[i]]
                xs = [data[key]['truthNeutronVtx_x_belle_frame'].iloc[i], data[key]['chipx'].iloc[i]]
                ax[1].plot(zs, xs, color = color[j], lw = 0.5, alpha = 0.5)
                j+=1
            ax[0].plot(data[key]['chipz'].mean(), data[key]['chipx'].mean(), 's', color = 'tab:red', markersize = 20, alpha = 0.5)
            ax[1].plot(data[key]['chipz'].mean(), data[key]['chipx'].mean(), 's', color = 'tab:red', markersize = 20, alpha = 0.5)
        
        plt.colorbar(sm, ax=ax[0]).set_label('Neutron Energy [keV]', rotation = 270, labelpad = 20)
        plt.colorbar(sm, ax=ax[1]).set_label('Neutron Energy [keV]', rotation = 270, labelpad = 20)
        #plt.colorbar(sm, ax=ax[0]).set_label(r'$E_{recoil}/E_{neutron}$', rotation = 270, labelpad = 20)
        #plt.colorbar(sm, ax=ax[1]).set_label(r'$E_{recoil}/E_{neutron}$', rotation = 270, labelpad = 20)
        
        plt.savefig("/home/jeef/Pictures/all_SimHits_%s.png"%(bgtype))
        #plt.savefig("neutron_production_points_ratio.png")
        plt.clf()

    def visualize_3D(self):
        cm = matplotlib.cm.viridis
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
        color = cm(norm(df['truthNeutronEnergy']))
        ###

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        j=0
        for key in data.keys():
            ax.plot(data[key]['chipz'].mean(), data[key]['chipx'].mean(), data[key]['chipy'].mean(),'s', color = 'tab:red', markersize = 20)
            for i in range(0,len(data[key])):
                zs = [data[key]['truthNeutronVtx_z_belle_frame'].iloc[i], data[key]['chipz'].iloc[i]]
                xs = [data[key]['truthNeutronVtx_x_belle_frame'].iloc[i], data[key]['chipx'].iloc[i]]
                ys = [data[key]['truthNeutronVtx_y_belle_frame'].iloc[i], data[key]['chipy'].iloc[i]]
                ax.plot(zs, xs, ys, color = color[j])
                j+=1
        #ax.scatter(df['truthNeutronVtx_z_belle_frame'], df['truthNeutronVtx_x_belle_frame'], df['truthNeutronVtx_y_belle_frame'], c = color)
        ax.set_xlabel('z')
        ax.set_ylabel('x')
        ax.set_zlabel('y')
        plt.colorbar(sm).set_label('Neutron Energy [keV]')
        plt.show()
        
a = analysis()
#data, truth = a.get_MC_data()
for bgtype in ['Coulomb_LER_dynamic', 'Coulomb_HER_dynamic', 'Brems_LER_dynamic', 'twoPhoton_Lumi', 'RBB_Lumi', 'Touschek_LER', 'Touschek_HER']:
    a.visualize_MC_with_geometry(bgtype)
#a.visualize_3D()
#df = a.get_simulated_rates()
