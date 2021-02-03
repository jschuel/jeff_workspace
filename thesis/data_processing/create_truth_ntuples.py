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

class create_truth_ntuples:

    def __init__(self):
        self.merge_raw_ntuples('truthB')
        self.merge_raw_ntuples('truth_e')
        self.merge_raw_ntuples('truth_mu')
        self.process_and_combine_merged_ntuples()
        self.clean()

    def merge_raw_ntuples(self, treename, filedir = '~/data/thesis_data/raw_analysis_ntuples/'):
        print("Merging ntuples with tree %s"%(treename))
        ch = ROOT.TChain(treename)
        ch.Add(filedir+'evtgen*.root')
        ch.Merge(filedir + 'all_'+treename+'.root')
        print("Merged!")

    def process_and_combine_merged_ntuples(self, infiledir = '~/data/thesis_data/raw_analysis_ntuples/', outfiledir = '~/data/thesis_data/truth_B_and_lepton_ntuples/', num_events_per_file = 1000):
        print("Processing and combining truth ntuples")
        try:
            Bs = rp.read_root(infiledir + 'all_truthB.root', key = 'truthB')
            es = rp.read_root(infiledir + 'all_truth_e.root', key = 'truth_e')
            mus = rp.read_root(infiledir + 'all_truth_mu.root', key = 'truth_mu')
        except:
            print('Ntuples not combined. Need to call "merge_raw_ntuples" first to create these three files')
            return 0
        #log true indices and original event numbers
        Bs['truth_index'] = Bs.index.to_numpy()
        Bs['original_event_number'] = Bs['__event__'].to_numpy()
        es['truth_index'] = es.index.to_numpy()
        es['original_event_number'] = es['__event__'].to_numpy()
        mus['truth_index'] = mus.index.to_numpy()
        mus['original_event_number'] = mus['__event__'].to_numpy()

        #label correct event numbers for analysis
        @jit(nopython=True)
        def update_event_numbers(events_list, num_events = num_events_per_file):
            print("Updating event numbers. This may take a while")
            for i in range(1,len(events_list)):
                while events_list[i]<events_list[i-1]:
                    events_list[i]+=num_events
            return events_list

        Bs['__event__'] = update_event_numbers(Bs['__event__'].to_numpy())
        es['__event__'] = update_event_numbers(es['__event__'].to_numpy())
        mus['__event__'] = update_event_numbers(mus['__event__'].to_numpy())

        leptons = es.append(mus)
        leptons = leptons.sort_values(by=['__event__'])
        leptons.index = [i for i in range(0,len(leptons))]
        print("Building ROOT files")
        Bs.to_root(outfiledir + "%s_truth_events.root"%(len(Bs)), key='truthB')
        leptons.to_root(outfiledir + "%s_truth_events.root"%(len(Bs)), key='truth_lepton', mode = 'a')

    def clean(self, filedir = '/home/jeff/data/thesis_data/raw_analysis_ntuples/'):
        print("Removing intermediate files")
        os.remove(filedir + 'all_truthB.root')
        os.remove(filedir + 'all_truth_e.root')
        os.remove(filedir + 'all_truth_mu.root')

create_truth_ntuples()        

        
        
