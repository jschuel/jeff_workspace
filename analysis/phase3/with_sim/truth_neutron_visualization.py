import uproot as ur
import ROOT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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

    def get_MC_data(self, bgtype):
        tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        #bgtype = ['Coulomb_HER_base', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_LER_dynamic', 'Brems_HER_base', 'Brems_LER_base', 'Brems_HER_dynamic', 'Brems_LER_dynamic', 'Touschek_HER_all', 'Touschek_LER_all']
        #bgtype = 'RBB_Lumi'
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

        cm = matplotlib.cm.plasma
        norm = matplotlib.colors.LogNorm(vmin = 1e0, vmax = 1e5)
        #norm = matplotlib.colors.Normalize(vmin = 0)
        sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)

        data = self.get_MC_data(bgtype)
        print(data)
        
        ###Make Consistent Color Scale###
        df = pd.DataFrame()
        for key in data.keys():
            #data[key] = data[key].loc[(data[key]['truthNeutronVtx_z_belle_frame']>-870) & #Line is for RBB hotspot
               #(data[key]['truthNeutronVtx_z_belle_frame']<-750) & (data[key]['truthNeutronVtx_x_belle_frame']>20) & 
               #(data[key]['truthNeutronVtx_x_belle_frame']<60) & (np.abs(data[key]['truthNeutronVtx_y_belle_frame'])<20)]
            df = df.append(data[key])
        df.index = [i for i in range(0,len(df))]
        color = cm(norm(df.dropna()['truthNeutronEnergy']))
        ###
        
        #img = plt.imread("/home/jeef/Pictures/farbeamline_update.png")
        img = plt.imread("/home/jeef/Pictures/farbeamline_nocolor.png")
        fig, ax = plt.subplots(2,1,figsize = (12,10))
        fig.suptitle('%s neutron SimHits and vertices'%(bgtype), fontsize=16)
        ax[0].set_xticks([])
        ax[0].set_yticks([])
        ax[1].set_xticks([])
        ax[1].set_yticks([])
        ax[0].set_xlabel('z [cm]')
        ax[0].set_ylabel('x [cm]')
        ax[1].set_xlabel('z [cm]')
        ax[1].set_ylabel('x [cm]')
        #ax[0].imshow(np.flipud(img), origin = 'lower', extent = [-3400,3190, -440,421], aspect = 'auto')
        #ax[1].imshow(np.flipud(img), origin = 'lower', extent = [-3400,3190, -440,421], aspect = 'auto')
        ax[0].imshow(np.flipud(img), origin = 'lower', extent = [-3333,3142, -438,414], aspect = 'auto')
        ax[1].imshow(np.flipud(img), origin = 'lower', extent = [-3333,3142, -438,414], aspect = 'auto')
        ax[0].add_patch(Rectangle((-1415.5,196), 31, 10, color = 'green', alpha = 0.5)) #elepaio
        ax[0].add_patch(Rectangle((-815.5,191), 31, 10, color = 'green', alpha = 0.5)) #tako
        ax[0].add_patch(Rectangle((-578,-189), 31, 10, color = 'green', alpha = 0.5)) #palila
        ax[0].add_patch(Rectangle((641,-199.4), 31, 10, color = 'green', alpha = 0.5)) #iiwi
        ax[0].add_patch(Rectangle((1385,172), 31, 10, color = 'green', alpha = 0.5)) #nene
        ax[0].add_patch(Rectangle((1585,170), 31, 10, color = 'green', alpha = 0.5)) #humu
        ax[1].add_patch(Rectangle((-1415.5,196), 31, 10, color = 'green', alpha = 0.5)) #elepaio
        ax[1].add_patch(Rectangle((-815.5,191), 31, 10, color = 'green', alpha = 0.5)) #tako
        ax[1].add_patch(Rectangle((-578,-189), 31, 10, color = 'green', alpha = 0.5)) #palila
        ax[1].add_patch(Rectangle((641,-199.4), 31, 10, color = 'green', alpha = 0.5)) #iiwi
        ax[1].add_patch(Rectangle((1385,172), 31, 10, color = 'green', alpha = 0.5)) #nene
        ax[1].add_patch(Rectangle((1585,170), 31, 10, color = 'green', alpha = 0.5)) #humu
        p = ax[0].scatter(df['truthNeutronVtx_z_belle_frame'], df['truthNeutronVtx_x_belle_frame'], c= df['truthNeutronEnergy'], cmap = 'plasma', s=0.4, norm = matplotlib.colors.LogNorm(vmin = 1e0, vmax = 1e5))
        ax[1].set_prop_cycle('color', color)
        ax[1].plot(df.dropna()[['truthNeutronVtx_z_belle_frame','chipz']].T, df.dropna()[['truthNeutronVtx_x_belle_frame','chipx']].T, lw=0.3, alpha = 0.15)
        plt.colorbar(p,ax=ax[0]).set_label(r'$E_{neutron}$', rotation = 270, labelpad = 20)
        plt.colorbar(sm,ax=ax[1]).set_label(r'E_{neutron}$', rotation = 270, labelpad = 20)
                
        plt.savefig("/home/jeef/Pictures/all_SimHits_%s.png"%(bgtype))
        #plt.savefig("neutron_production_points_ratio.png")
        plt.clf()

    def visualize_3D(self, bgtype):
        cm = matplotlib.cm.viridis
        norm = matplotlib.colors.LogNorm(vmin = 1e1, vmax = 1e5)
        #norm = matplotlib.colors.Normalize(vmin = 0)
        sm = matplotlib.cm.ScalarMappable(cmap=cm, norm=norm)

        data = self.get_MC_data(bgtype)

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
            ax.plot(data[key]['chipz'].mean(), data[key]['chipx'].mean(), data[key]['chipy'].mean(),'s', color = 'tab:red', markersize = 20, alpha = 0.5)
            for i in range(0,len(data[key])):
                zs = [data[key]['truthNeutronVtx_z_belle_frame'].iloc[i], data[key]['chipz'].iloc[i]]
                xs = [data[key]['truthNeutronVtx_x_belle_frame'].iloc[i], data[key]['chipx'].iloc[i]]
                ys = [data[key]['truthNeutronVtx_y_belle_frame'].iloc[i], data[key]['chipy'].iloc[i]]
                ax.plot(zs, xs, ys, 'o', color = color[j], markersize = 0.8)
                j+=1
        #ax.scatter(df['truthNeutronVtx_z_belle_frame'], df['truthNeutronVtx_x_belle_frame'], df['truthNeutronVtx_y_belle_frame'], c = color)
        ax.set_xlabel('z')
        ax.set_ylabel('x')
        ax.set_zlabel('y')
        plt.colorbar(sm).set_label('Neutron Energy [keV]')
        fig.suptitle('%s 3D neutron SimHits and vertices'%(bgtype), fontsize=16)
        plt.savefig("/home/jeef/Pictures/3D_all_SimHits_%s.png"%(bgtype))
        plt.show()

    def get_all_MC_data(self):
        tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        bgtype = ['Coulomb_LER_dynamic', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_HER_base', 'Brems_LER_dynamic', 'Brems_HER_dynamic', 'Brems_LER_base', 'twoPhoton_Lumi', 'RBB_Lumi', 'Touschek_LER', 'Touschek_HER']
        truth = {}
        for tpc in tpcs:
            dir = '/home/jeef/data/phase3/spring_2020/05-09-20/geant4_simulation/truth_neutrons_only/%s/'%(tpc)
            for bg in bgtype:
                try:
                    try:
                        truth[tpc+'_'+bg] = ur.open(dir+bg+'_'+tpc+'_all.root')['recoils'].pandas.df(flatten=False)
                    except KeyError:
                        truth[tpc+'_'+bg] = ur.open(dir+bg+'_'+tpc+'_all.root')['my_ttree'].pandas.df(flatten=False)
                    try:
                        truth[tpc+'_'+bg]['chipx'] = truth[tpc+'_'+bg]['chipx'].apply(lambda x: x[0])
                        truth[tpc+'_'+bg]['chipy'] = truth[tpc+'_'+bg]['chipy'].apply(lambda x: x[0])
                        truth[tpc+'_'+bg]['chipz'] = truth[tpc+'_'+bg]['chipz'].apply(lambda x: x[0])
                    except TypeError:
                        pass
                except FileNotFoundError:
                    truth[tpc+'_'+bg] = pd.DataFrame()
        return truth

    def visualize_all_with_geometry(self):
        plt.rc('legend', fontsize=12)
        plt.rc('xtick', labelsize=16)
        plt.rc('ytick', labelsize=16)
        plt.rc('axes', labelsize=18)
        plt.rc('axes', titlesize=18)

        bgtypes = ['Coulomb_LER_dynamic', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_HER_base', 'Brems_LER_dynamic', 'Brems_HER_dynamic', 'Brems_LER_base', 'twoPhoton_Lumi', 'RBB_Lumi', 'Touschek_LER', 'Touschek_HER']
        data = self.get_all_MC_data()
        df = pd.DataFrame()
        for key in data.keys():
            df = df.append(data[key])
        df.index = [i for i in range(0,len(df))]
        df_iiwi = df.loc[df['detNb'] == 2]
        df_palila = df.loc[df['detNb'] == 3]
        df_tako = df.loc[df['detNb'] == 4]
        df_elepaio = df.loc[df['detNb'] == 5]
        #img = plt.imread("/home/jeef/Pictures/farbeamline_update.png")
        img = plt.imread("/home/jeef/Pictures/farbeamline_nocolor.png")
        fig, ax = plt.subplots(1,1,figsize = (25,10))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('z [cm]')
        ax.set_ylabel('x [cm]')
        #ax.imshow(np.flipud(img), origin = 'lower', extent = [-3400,3190, -440,421], aspect = 'auto')
        ax.imshow(np.flipud(img), origin = 'lower', extent = [-3333,3142, -438,414], aspect = 'auto')
        #p = ax.scatter(df['truthNeutronVtx_z_belle_frame'],df['truthNeutronVtx_x_belle_frame'],c=df['truthNeutronEnergy'],s=0.8,norm = matplotlib.colors.LogNorm(vmin = 1e1, vmax = 1e5),cmap = 'plasma')
        #ax.plot(df_iiwi['chipz'].mean(), df_iiwi['chipx'].mean() ,'s', color = 'green', markersize = 8, alpha = 0.5)
        #ax.plot(df_palila['chipz'].mean(), df_palila['chipx'].mean(), 's', color = 'green', markersize = 8, alpha = 0.5)
        #ax.plot(df_tako['chipz'].mean(), df_tako['chipx'].mean(),'s', color = 'green', markersize = 8, alpha = 0.5)
        #ax.plot(df_elepaio['chipz'].mean(), df_elepaio['chipx'].mean()-5,'s', color = 'green', markersize = 8, alpha = 0.5)

        ax.add_patch(Rectangle((-1415.5,196), 31, 10, color = 'green', alpha = 0.5)) #elepaio
        ax.add_patch(Rectangle((-815.5,191), 31, 10, color = 'green', alpha = 0.5)) #tako
        ax.add_patch(Rectangle((-578,-189), 31, 10, color = 'green', alpha = 0.5)) #palila
        ax.add_patch(Rectangle((641,-199.4), 31, 10, color = 'green', alpha = 0.5)) #iiwi
        ax.add_patch(Rectangle((1385,172), 31, 10, color = 'green', alpha = 0.5)) #nene
        ax.add_patch(Rectangle((1585,170), 31, 10, color = 'green', alpha = 0.5)) #humu

        ### Old
        #ax.plot([-1400,-800,-562.5], [204, 204,-181.4],'s', color = 'green', markersize = 6, alpha = 0.5)
        #ax.plot([656.5, 1400, 1600], [-182.9, 201,201], 's', color = 'green', markersize = 6, alpha = 0.5)
        #ax.plot([-1400,-800,-562.5], [201,196,-184.0],'s', color = 'green', markersize = 6, alpha = 0.5)
        #ax.plot([656.5, 1400, 1600], [-194.4, 177,175], 's', color = 'green', markersize = 6, alpha = 0.5)
        ### End old
        
        #plt.colorbar(p).set_label('Neutron Energy [keV]')
        #fig.suptitle('All neutron SimHits and vertices', fontsize=16)
        #plt.savefig("/home/jeef/Pictures/3D_all_SimHits_%s.png"%(bgtype))
        #plt.savefig("/home/jeef/Pictures/old.png")
        plt.savefig("/home/jeef/Pictures/new.png")
        plt.show()

    def visualize_all_3D(self):
        bgtypes = ['Coulomb_LER_dynamic', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_HER_base', 'Brems_LER_dynamic', 'Brems_HER_dynamic', 'Brems_LER_base', 'twoPhoton_Lumi', 'RBB_Lumi', 'Touschek_LER', 'Touschek_HER']
        data = self.get_all_MC_data()
        df = pd.DataFrame()
        for key in data.keys():
            df = df.append(data[key])
        df.index = [i for i in range(0,len(df))]
        df_iiwi = df.loc[df['detNb'] == 2]
        df_palila = df.loc[df['detNb'] == 3]
        df_tako = df.loc[df['detNb'] == 4]
        df_elepaio = df.loc[df['detNb'] == 5]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        #j=0
        
        #for i in range(0,len(df)):
        #    zs = [df['truthNeutronVtx_z_belle_frame'].iloc[i], df['chipz'].iloc[i]]
        #    xs = [df['truthNeutronVtx_x_belle_frame'].iloc[i], df['chipx'].iloc[i]]
        #    ys = [df['truthNeutronVtx_y_belle_frame'].iloc[i], df['chipy'].iloc[i]]
        #    ax.plot(zs, xs, ys, 'o', color = color[j], markersize = 0.8)
        #    j+=1
        p = ax.scatter(df['truthNeutronVtx_z_belle_frame'],df['truthNeutronVtx_x_belle_frame'],df['truthNeutronVtx_y_belle_frame'],c=df['truthNeutronEnergy'],s=0.8,norm = matplotlib.colors.LogNorm(vmin = 1e1, vmax = 1e5), cmap = 'plasma')
        ax.plot(df_iiwi['chipz'].mean(), df_iiwi['chipx'].mean(), df_iiwi['chipy'].mean(),'s', color = 'green', markersize = 20, alpha = 0.5)
        ax.plot(df_palila['chipz'].mean(), df_palila['chipx'].mean(), df_palila['chipy'].mean(),'s', color = 'green', markersize = 20, alpha = 0.5)
        ax.plot(df_tako['chipz'].mean(), df_tako['chipx'].mean(), df_tako['chipy'].mean(),'s', color = 'green', markersize = 20, alpha = 0.5)
        ax.plot(df_elepaio['chipz'].mean(), df_elepaio['chipx'].mean(), df_elepaio['chipy'].mean(),'s', color = 'green', markersize = 20, alpha = 0.5)
        ax.plot([1400, 1600], [201,201], [16, 16], 's', color = 'green', markersize = 20, alpha = 0.5)
        ax.set_xlabel('z')
        ax.set_ylabel('x')
        ax.set_zlabel('y')
        plt.colorbar(p).set_label('Neutron Energy [keV]')
        fig.suptitle('All neutron SimHits and vertices', fontsize=16)
        #plt.savefig("/home/jeef/Pictures/3D_all_SimHits_%s.png"%(bgtype))
        plt.show()
        
    
a = analysis()
#data, truth = a.get_MC_data()
for bgtype in ['Touschek_LER']: #['Coulomb_LER_dynamic', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_HER_base', 'Brems_LER_dynamic', 'Brems_HER_dynamic', 'Brems_LER_base', 'twoPhoton_Lumi', 'RBB_Lumi', 'Touschek_LER', 'Touschek_HER']:
    a.visualize_MC_with_geometry(bgtype)
    #a.visualize_3D(bgtype)
#df = a.get_simulated_rates()
#a.visualize_all_3D()
#a.visualize_all_with_geometry()
