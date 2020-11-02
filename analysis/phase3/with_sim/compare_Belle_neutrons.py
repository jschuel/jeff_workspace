import uproot as ur
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class compare_neutron_rates:

    def __init__(self):
        self.nfbl, self.fbl = self.load_files()

    def load_files(self, base_path = '/home/jeef/data/phase3/Belle2_neutron_ntuples/'):
        farbeamline_path = base_path + 'farbeamline/'
        nofarbeamline_path = base_path + 'no_farbeamline/'
        dfs_fbl = {}
        dfs = {}
        for f in os.listdir(farbeamline_path):
            dfs_fbl[f] = ur.open(farbeamline_path + f)['tree'].pandas.df()
        for f in os.listdir(nofarbeamline_path):
            dfs[f] = ur.open(nofarbeamline_path + f)['tree'].pandas.df()

        return dfs, dfs_fbl

    def compare_fractions_of_events(self, detector):
        detector = detector.upper()
        detID = {'IR':0, 'PXD':1, 'SVD':2, 'CDC':3, 'ARICH':4, 'TOP':5, 'ECL':6, 'EKLM':7, 'BKLM':8}
        fbl, nfbl = self.fbl, self.nfbl
        for key in fbl.keys():
            neutrons_inside = np.abs(nfbl[key].loc[nfbl[key]['subDet'] == detID[detector]]['neutronWeight']).to_numpy().sum()
            neutrons_outside = np.abs(fbl[key].loc[fbl[key]['subDet'] == detID[detector]]['neutronWeight']).to_numpy().sum() - neutrons_inside
            frac = neutrons_outside/(neutrons_inside+neutrons_outside)
            print(f, frac)

c = compare_neutron_rates()
c.compare_fractions_of_events('BKLM')
