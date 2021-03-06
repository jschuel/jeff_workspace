''' JTS 02/2/21
Starts with raw MC analysis ntuples that were created with basf2's b2bii module and
turns them into processed combined analysis ntuples that include truth B and truth lepton
information. The general procedure is as follows:

1. Use TChain to concatenate all ntuples of an individual branch together. Do this for
truthB, truth_e, and truth_mu branches
2. Save these combined files as intermediate files
3. Load intermediate files as dataframes
4. Save truth indices and relabel event numbers accordingly (since they were created in batches)
5. Merge dataframes into a root file containing truthB and truth_lepton branches
6. Clean up and delete the intermediate files
'''

import numpy as np
import pandas as pd
import ROOT
import root_pandas as rp
import os
from numba import jit

class merge_ntuples:

    def __init__(self):
        #for key in ['B', 'Dst', 'D', 'K', 'pi0', 'pi', 'e', 'mu']:
        for key in ['B_truth', 'e_truth', 'mu_truth', 'e','mu']:
            self.merge_raw_ntuples(key)
        self.process_and_combine_merged_ntuples()
        #self.clean()

    def merge_raw_ntuples(self, treename, filedir = '~/data/thesis_data/raw_reco_ntuples/'):
        print("Merging ntuples with tree %s"%(treename))
        ch = ROOT.TChain(treename)
        ch.Add(filedir+'evtgen*.root')
        ch.Merge(filedir + 'all_'+treename+'.root')
        print("Merged!")

    def process_and_combine_merged_ntuples(self, infiledir = '~/data/thesis_data/raw_reco_ntuples/', outfiledir = '~/data/thesis_data/combined_reco_ntuples/', num_events_per_file = 10000, keys = ['B_truth', 'e_truth', 'mu_truth', 'e','mu']):
        print("Processing and combining truth ntuples")
        try:
            dfs = {}
            for key in keys:
                dfs[key] = rp.read_root(infiledir + 'all_%s.root'%(key), key = key)
                print("Read in all_%s.root"%(key))
        except:
            print('Error: Check file names')
            return 0
        #log true indices and original event numbers
        for key in keys:
            dfs[key]['truth_index'] = dfs[key].index.to_numpy()
            dfs[key]['original_event_number'] = dfs[key]['__event__'].to_numpy()

        #label correct event numbers for analysis
        @jit(nopython=True)
        def update_event_numbers(events_list, num_events = num_events_per_file):
            print("Updating event numbers. This may take a while")
            for i in range(1,len(events_list)):
                while events_list[i]<events_list[i-1]:
                    events_list[i]+=num_events
            return events_list

        
        for i, key in enumerate(keys):
            dfs[key]['__event__'] = update_event_numbers(dfs[key]['__event__'].to_numpy())
            dfs[key] = dfs[key].sort_values(by=['__event__'])
            dfs[key].index = [i for i in range(0,len(dfs[key]))]
            print("Building ROOT files")
            if i == 0:
                dfs[key].to_root(outfiledir + "merged_leptons2.root", key=key)
            else:
                dfs[key].to_root(outfiledir + "merged_leptons2.root", key=key, mode='a')

    def clean(self, filedir = '/home/jeff/data/thesis_data/raw_analysis_ntuples/'):
        print("Removing intermediate files")
        #os.remove(filedir + 'all_e.root')
        #os.remove(filedir + 'all_mu.root')

merge_ntuples()        

        
        
