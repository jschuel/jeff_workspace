import root_pandas as rp
import uproot as ur
import ROOT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

class extract:

    def __init__(self):
        self.tpc_list = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        pass

    def get_MC_data(self, bgtype='RBB_Lumi'):
        tpcs = self.tpc_list
        truth = {}
        for tpc in tpcs:
            dir = '~/data/phase3/spring_2020/05-09-20/geant4_simulation/truth_neutrons_only/%s/'%(tpc)
            try:
                try:
                    truth[tpc] = rp.read_root(dir+bgtype+'_'+tpc+'_all.root')
                except KeyError:
                    truth[tpc] = rp.read_root(dir+bgtype+'_'+tpc+'_all.root')
                try:
                    truth[tpc]['chipx'] = truth[tpc]['chipx'].apply(lambda x: x[0])
                    truth[tpc]['chipy'] = truth[tpc]['chipy'].apply(lambda x: x[0])
                    truth[tpc]['chipz'] = truth[tpc]['chipz'].apply(lambda x: x[0])
                except TypeError:
                    pass
            except FileNotFoundError:
                truth[tpc] = pd.DataFrame()
        return truth

    def apply_RBB_cuts(self, bgtype='RBB_Lumi'):
        data = self.get_MC_data(bgtype=bgtype)
        tpcs = self.tpc_list
        dfs = {}
        for tpc in tpcs:
            if tpc == 'iiwi' or tpc == 'nene' or tpc =='humu':
                #dfs[tpc] = data[tpc].loc[(data[tpc]['truthNeutronVtx_z_belle_frame']>1350) & (data[tpc]['truthNeutronVtx_z_belle_frame']<1900) 
                #& (data[tpc]['truthNeutronVtx_x_belle_frame']<95) & (data[tpc]['truthNeutronVtx_x_belle_frame']>37) 
                #& (np.abs(data[tpc]['truthNeutronVtx_y_belle_frame'])<20)]
                dfs[tpc] = data[tpc].loc[(data[tpc]['truthNeutronVtx_z_belle_frame']>=1350) & (data[tpc]['truthNeutronVtx_z_belle_frame']<=1650) 
                & (data[tpc]['truthNeutronVtx_x_belle_frame']<=80) & (data[tpc]['truthNeutronVtx_x_belle_frame']>=40) 
                & (np.abs(data[tpc]['truthNeutronVtx_y_belle_frame'])<=20)]
            else:
                dfs[tpc] = data[tpc].loc[(data[tpc]['truthNeutronVtx_z_belle_frame']>=-900) & 
               (data[tpc]['truthNeutronVtx_z_belle_frame']<=-750) & (data[tpc]['truthNeutronVtx_x_belle_frame']>=16) & 
               (data[tpc]['truthNeutronVtx_x_belle_frame']<=60) & (np.abs(data[tpc]['truthNeutronVtx_y_belle_frame'])<=20)]
            dfs[tpc].index = [i for i in range(0, len(dfs[tpc]))]
        return dfs

    ### Code below is for making paper figure with all vertices ###

    def get_Lumi_MC_data(self):
        tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        bgtype = ['twoPhoton_Lumi', 'RBB_Lumi']
        truth = {}
        for tpc in tpcs:
            dir = '~/data/phase3/spring_2020/05-09-20/geant4_simulation/truth_neutrons_only/%s/'%(tpc)
            for bg in bgtype:
                try:
                    try:
                        #truth[tpc+'_'+bg] = ur.open(dir+bg+'_'+tpc+'_all.root')['recoils'].pandas.df(flatten=False)
                        truth[tpc+'_'+bg] = rp.read_root(dir+bg+'_'+tpc+'_all.root', key='recoils')
                    except KeyError:
                        #truth[tpc+'_'+bg] = ur.open(dir+bg+'_'+tpc+'_all.root')['my_ttree'].pandas.df(flatten=False)
                        truth[tpc+'_'+bg] = rp.read_root(dir+bg+'_'+tpc+'_all.root', key='my_ttree')
                    try:
                        truth[tpc+'_'+bg]['chipx'] = truth[tpc+'_'+bg]['chipx'].apply(lambda x: x[0])
                        truth[tpc+'_'+bg]['chipy'] = truth[tpc+'_'+bg]['chipy'].apply(lambda x: x[0])
                        truth[tpc+'_'+bg]['chipz'] = truth[tpc+'_'+bg]['chipz'].apply(lambda x: x[0])
                    except TypeError:
                        pass
                except FileNotFoundError:
                    truth[tpc+'_'+bg] = pd.DataFrame()
        return truth

    def get_single_beam_MC_data(self):
        tpcs = ['iiwi', 'nene', 'humu', 'palila', 'tako', 'elepaio']
        bgtype = ['Coulomb_LER_dynamic', 'Coulomb_LER_base', 'Coulomb_HER_dynamic', 'Coulomb_HER_base', 'Brems_LER_dynamic', 'Brems_HER_dynamic', 'Brems_LER_base', 'Touschek_LER', 'Touschek_HER']
        truth = {}
        for tpc in tpcs:
            dir = '~/data/phase3/spring_2020/05-09-20/geant4_simulation/truth_neutrons_only/%s/'%(tpc)
            for bg in bgtype:
                try:
                    try:
                        #truth[tpc+'_'+bg] = ur.open(dir+bg+'_'+tpc+'_all.root')['recoils'].pandas.df(flatten=False)
                        truth[tpc+'_'+bg] = rp.read_root(dir+bg+'_'+tpc+'_all.root',key='recoils')
                    except KeyError:
                        #truth[tpc+'_'+bg] = ur.open(dir+bg+'_'+tpc+'_all.root')['my_ttree'].pandas.df(flatten=False)
                        truth[tpc+'_'+bg] = rp.read_root(dir+bg+'_'+tpc+'_all.root',key='my_ttree')

                    try:
                        truth[tpc+'_'+bg]['chipx'] = truth[tpc+'_'+bg]['chipx'].apply(lambda x: x[0])
                        truth[tpc+'_'+bg]['chipy'] = truth[tpc+'_'+bg]['chipy'].apply(lambda x: x[0])
                        truth[tpc+'_'+bg]['chipz'] = truth[tpc+'_'+bg]['chipz'].apply(lambda x: x[0])
                    except TypeError:
                        pass
                except FileNotFoundError:
                    truth[tpc+'_'+bg] = pd.DataFrame()
        return truth

    def make_paper_figure(self):
        plt.rc('legend', fontsize=20)
        plt.rc('xtick', labelsize=24)
        plt.rc('ytick', labelsize=24)
        plt.rc('axes', labelsize=26)
        plt.rc('axes', titlesize=26)
        
        lumi = self.get_Lumi_MC_data()
        sb = self.get_single_beam_MC_data()
        
        df_SB = pd.DataFrame()
        df_Lumi = pd.DataFrame()
        for key in lumi.keys():
            df_Lumi = df_Lumi.append(lumi[key])
        for key in sb.keys():
            df_SB = df_SB.append(sb[key])

        df_SB = df_SB.loc[np.abs(df_SB['truthNeutronVtx_z_belle_frame'])>5]
        df_Lumi = df_Lumi.loc[np.abs(df_Lumi['truthNeutronVtx_z_belle_frame'])>5]
        
        df_SB.index = [i for i in range(0,len(df_SB))]
        df_Lumi.index = [i for i in range(0,len(df_Lumi))]
        img = plt.imread("/home/jeff/Pictures/farbeamline_nocolor.png")
        fig, ax = plt.subplots(1,1,figsize = (25,10))
        ax.set_xlabel('z [cm]')
        ax.set_ylabel('x [cm]')
        ax.set_ylim(-300,300)
        ax.set_xlim(-2800,2800)
        ax.imshow(np.flipud(img), origin = 'lower', extent = [-3333,3142, -438,414], aspect = 'auto')

        ax.plot(df_SB[['truthNeutronVtx_z_belle_frame','chipz']].T, df_SB[['truthNeutronVtx_x_belle_frame','chipx']].T, lw=0.15, alpha = 0.15, color = 'magenta')

        ax.plot(df_Lumi[['truthNeutronVtx_z_belle_frame','chipz']].T, df_Lumi[['truthNeutronVtx_x_belle_frame','chipx']].T, lw=0.15, alpha = 0.15, color = 'cyan')

        ax.add_patch(Rectangle((-1415.5,196), 31, 10, facecolor = 'gold', alpha = 1.0, edgecolor = 'black', zorder=1e6)) #elepaio
        ax.add_patch(Rectangle((-815.5,191), 31, 10, facecolor = 'gold', alpha = 1.0, edgecolor = 'black', zorder=1e6+1)) #tako
        ax.add_patch(Rectangle((-578,-189), 31, 10, facecolor = 'gold', alpha = 1.0, edgecolor = 'black', zorder=1e6+2)) #palila
        ax.add_patch(Rectangle((641,-199.4), 31, 10, facecolor = 'gold', alpha = 1.0, edgecolor = 'black', zorder=1e6+3)) #iiwi
        ax.add_patch(Rectangle((1385,172), 31, 10, facecolor = 'gold', alpha = 1.0, edgecolor = 'black', zorder=1e6+4)) #nene
        ax.add_patch(Rectangle((1585,170), 31, 10, facecolor = 'gold', alpha = 1.0, edgecolor = 'black', zorder=1e6+5)) #humu

        ax.add_patch(Rectangle((-900,18), 150, 40, facecolor = 'lime', alpha = 0.7, edgecolor = 'black', zorder=1e6+6)) #hotspot
        ax.add_patch(Rectangle((1350,40), 300, 40, facecolor = 'lime', alpha = 0.7, edgecolor = 'black', zorder=1e6+7)) #hotspot

        handle = [Line2D([0], [0], marker='s', color='black', label='TPC',
                          markerfacecolor='gold', alpha = 1.0, markersize=15, lw=0),
        Line2D([0], [0], color='magenta', label='Single Beam', alpha = 1.0),
        Line2D([0], [0], color = 'cyan', label='Luminosity', alpha = 1.0)]
    
        plt.legend(handles=handle, bbox_to_anchor=(0.84, 0.8), framealpha=1)
        plt.savefig("/home/jeff/Pictures/MC_tracks_update.png")
        plt.show()
        
a = extract()
