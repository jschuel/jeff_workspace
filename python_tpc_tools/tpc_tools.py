import root_pandas as rp
import pandas as pd
import numpy as np
import array
import matplotlib.pyplot as plt
import sklearn.cluster as cluster
from ROOT import TVector3
import ROOT
import re

class tpc_tools:
    def __init__(self, input_file='4500_honu_stop_mode_ext_trigger_scan_interpreted.h5', module_id="honu", output_file = '4500_honu.root'):
        self.calibration = self.get_calibration_info(module_id)
        self.hits = self.get_hits_data(input_file, module_id)
        self.meta = self.get_meta_data(input_file)
        self.data = self.merge_hits_and_meta_data(input_file, module_id)
        self.tracks = self.process_tracks()
        self.make_ntuple(output_file, alpha=10)

    def get_calibration_info(self, module_id): #Consider updating to yaml headers
        df = pd.read_pickle("/Users/vahsengrouplaptop/workspace/duke_testbeam/tpc_tools/1_initial_processing/calibration_input_files/calibrated_files/calibrated_post_phase2_%s.pkl"%(module_id))
        return df

    def get_hits_data(self, input_file, module_id):
        calibration = self.get_calibration_info(module_id)
        df_hits = pd.read_hdf(input_file, "Hits")
        df_hits = df_hits.drop(columns= ['trigger_number', 'trigger_time_stamp', 'TDC', 'TDC_time_stamp', 'trigger_status', 'service_record', 'event_status', 'LVL1ID', 'BCID'])
        df_hits = df_hits.rename(columns = {'relative_BCID': 'BCID'})
        df_hits['pixel_charge']=df_hits['tot'].map(calibration.set_index('tot_code')['q_per_tot'])
        df_hits['pixel_energy']=df_hits['tot'].map(calibration.set_index('tot_code')['E_per_tot'])
        return df_hits
         
    def get_meta_data(self, input_file):
        df_meta = pd.read_hdf(input_file, "meta_data") #generates dataframe of metadata
        df_meta = df_meta.drop(columns = ['error_code'])
        return df_meta

    def merge_hits_and_meta_data(self, input_file, module_id):
        hits = self.get_hits_data(input_file, module_id)
        meta = self.get_meta_data(input_file)
        meta = meta.loc[meta['event_number'].isin(hits['event_number'])==True] #find events numbers that match
        meta = meta.loc[meta.duplicated(['event_number']) == False] #remove duplicate numbers
        meta.index = [i for i in range(0,len(meta))]
        hits = hits.loc[hits['event_number'].isin(meta['event_number'])==True] #match up hits numbers to data numbers
        hits.index = [i for i in range(0,len(hits))]
        hits['column'], hits['row'] = hits['column']-1, hits['row']-1
        merged_data = hits.groupby('event_number')[[col for col in hits.columns if col != 'event_number']].agg(list) #creates one entry per event number where hit level data is stored as lists
        merged_data['raw_event_number'] = merged_data.index.to_numpy()
        merged_data.index = [i for i in range(0,len(merged_data))]
        merged_data['timestamp_start'] = meta['timestamp_start']
        merged_data['timestamp_stop'] = meta['timestamp_stop']
        merged_data['event_number'] = merged_data.index.to_numpy()
        merged_data['nhits'] = merged_data['column'].apply(lambda x: len(x))
        merged_data = merged_data[['event_number', 'raw_event_number', 'timestamp_start', 'timestamp_stop', 'nhits', 'column', 'row', 'tot', 'BCID', 'pixel_charge', 'pixel_energy']]
        for col in ['column', 'row', 'tot', 'BCID', 'pixel_charge', 'pixel_energy']:
            merged_data[col] = merged_data[col].apply(lambda x: np.array(x)) #turns lists into numpy arrays
        merged_data['sum_tot'] = merged_data['tot'].apply(lambda x: x.sum()) #add new columns
        merged_data['track_charge'] = merged_data['pixel_charge'].apply(lambda x: x.sum())
        merged_data['track_energy'] = merged_data['pixel_energy'].apply(lambda x: x.sum())
        merged_data['x'], merged_data['y'], merged_data['z'] = self.get_xyz_from_col_row_bcid(merged_data['column'], merged_data['row'], merged_data['BCID'])
        return merged_data

    def process_tracks(self):
        data = self.data
        data = self.get_hitside(data)
        data = data.loc[data['nhits']<10000]
        data.index = [i for i in range(0,len(data))]
        length = []
        theta = []
        phi = []
        head_q = []
        tail_q = []
        for i in range(0,len(data)):
            
            l, t, p, hq, tq = self.fit_track(data['x'][i], data['y'][i], data['z'][i])
            length.append(l)
            theta.append(t)
            phi.append(p)
            head_q.append(hq)
            tail_q.append(tq)
            print('fit track %s out of %s'%(i, len(data)))
        data['length'], data['theta'], data['phi'], data['head_charge'], data['tail_charge'] = length, theta, phi, head_q, tail_q
        neutron_cut = "track_energy < (0.5*length-75) and track_energy > (0.04*length-65) and track_energy > 100 and hitside_col_min == 0 and hitside_col_max == 0 and hitside_row_min == 0 and hitside_row_max == 0"
        alpha_cut = "hitside_col_min == 1 and hitside_col_max == 1 and hitside_row_min == 0 and hitside_row_max == 0 and track_energy > 800"
        xray_cut = "hitside_col_min == 0 and hitside_col_max == 0 and hitside_row_min == 0 and hitside_row_max == 0 and track_energy < 50 and length > 3000"
        data['is_tight_neutron'] = 0
        data['is_tight_alpha'] = 0
        data['is_tight_xray'] = 0
        data['is_not_identified'] = 0
        data['is_tight_neutron'][data.query(neutron_cut).index.to_numpy()]=1
        data['is_tight_alpha'][data.query(alpha_cut).index.to_numpy()]=1
        data['is_tight_xray'][data.query(xray_cut).index.to_numpy()]=1
        data['is_not_identified'][data.loc[(data['is_tight_neutron'] == 0) & (data['is_tight_alpha'] == 0) & (data['is_tight_xray'] == 0)].index.to_numpy()] = 1
        
        return data

    def make_quick_ntuple(self, output_file):
        f = re.search(r'(.*).root', output_file).group(1) #finds the non .root part of file name
        filename = f + '_tracks_only.root'
        data = self.tracks
        data.to_root(filename, key = "data") #use pandas to_root to process track level quantities

    def make_ntuple(self, output_file, **kwargs):
        df = self.tracks
        f = re.search(r'(.*).root', output_file).group(1) #finds the non .root part of file name
        key = [key for key in kwargs.keys()][0]
        filename = f + '_%s_%s_events.root'%(key, kwargs[key])
        if len(kwargs) == 1:
            df = df.loc[df['is_tight_%s'%(key)] == 1]
            if kwargs[key] >= len(df):
                print("There are only %s %s's in the file. Processesing all %s's"%(kwargs[key], key, key))
            elif kwargs[key] < len(df):
                df = df.head(kwargs[key])
        df.index = [i for i in range(0,len(df))]
        keys = [val for val in df.columns]
        output = ROOT.TFile(filename, 'recreate')
        tout = ROOT.TTree('data','data')
        branches = {}
        data={}
        
        for key in keys:
            if df[key].dtype == "O": #Determines the size of an array in a dataframe to be pushed to ntuple
                data[key]=array.array('d',[0 for j in range(0,50000)])
                branches[key]=tout.Branch("%s"%(key), data[key], "%s[nhits]/D"%(key))
            elif key == 'nhits':
                data[key]=array.array('i',[0])
                branches[key]=tout.Branch("%s"%(key), data[key], "%s/I"%(key))
            else:
                data[key]=array.array('d',[0])
                branches[key]=tout.Branch("%s"%(key), data[key], "%s/D"%(key))

        for j in range(0,len(df)):
            data['nhits'][0] = df['nhits'].to_numpy()[j].astype(int)
            for key in keys:
                if df[key].dtype == "O":
                    for i in range(0,data['nhits'][0]):
                        data[key][i]=df[key][j][i]
                elif key != 'nhits':
                    data[key][0]=df[key][j]
            print("Filling event %s out of %s"%(j, len(df)))
            tout.Fill()

        output.Write()
        output.Close()
    
    def get_xyz_from_col_row_bcid(self, col, row, bcid):
        return col*250., (335-row)*50., 250.*bcid

    def get_hitside(self, df): #Computes the edgecut information
        df['col_min'], df['col_max'], df['row_min'], df['row_max'], df['BCID_max'] = [df['column'][j].min() for j in range(0,len(df))], [df['column'][j].max() for j in range(0,len(df))], [df['row'][j].min() for j in range(0,len(df))], [df['row'][j].max() for j in range(0,len(df))], [df['BCID'][j].max() for j in range(0,len(df))]
        df['hitside_col_min'], df['hitside_col_max'], df['hitside_row_min'], df['hitside_row_max'], df['hitside_BCID_max'] = [0 for i in range(0,len(df))], [0 for i in range(0,len(df))], [0 for i in range(0,len(df))], [0 for i in range(0,len(df))], [0 for i in range(0,len(df))]
        df['hitside_col_min'][df.loc[df['col_min'] == 0].index.to_numpy()] = 1
        df['hitside_col_max'][df.loc[df['col_max'] == 79].index.to_numpy()] = 1
        df['hitside_row_min'][df.loc[df['row_min'] == 0].index.to_numpy()] = 1
        df['hitside_row_max'][df.loc[df['row_max'] == 335].index.to_numpy()] = 1
        df['hitside_BCID_max'][df.loc[df['BCID_max'] == 99].index.to_numpy()] = 1
        df = df.drop(columns = ['col_min', 'col_max', 'row_min', 'row_max', 'BCID_max'])
        return df

    def fit_track(self, x,y,z): #faster than root fit, so this is standard
        data = np.concatenate((x[:, np.newaxis], y[:, np.newaxis], z[:, np.newaxis]), axis=1)
        datamean = data.mean(axis = 0)
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
        head_q = 0
        tail_q = 0
        for p in projection:
            if p > midp:
                head_charge += 1
                head_q += 1
            else:
                tail_charge += 1
                tail_q += 1
        if head_charge > tail_charge:
            vec = -1.*vec
            tail_charge = head_q
            head_charge = tail_q
        head_charge_fraction = head_charge/(head_charge+tail_charge)
        tail_charge_fraction = 1 - head_charge_fraction
        theta = vec.Theta()*180/np.pi
        phi = vec.Phi()*180/np.pi
        return length, theta, phi, head_charge_fraction, tail_charge_fraction
        

        #return hits
        
class tpc_track(tpc_tools):
    def __init__(self, input_file='4500_honu_stop_mode_ext_trigger_scan_interpreted.h5', module_id="honu"):
        super().__init__(input_file='4500_honu_stop_mode_ext_trigger_scan_interpreted.h5', module_id="honu")

    def process_track(self,event_number):
        track = self.data.loc[self.data['event_number']==event_number]
        track['length'], track['theta'], track['phi'], track['head_charge'], track['tail_charge'] = self.fit_track(track['x'].iloc[0], track['y'].iloc[0], track['z'].iloc[0])
        
        ###hitside info###
        
        if track['column'].iloc[0].min() == 0:
            track['hitside_col_min'] = 1
        else:
            track['hitside_col_min'] = 0
        if track['column'].iloc[0].max() == 79:
            track['hitside_col_max'] = 1
        else:
            track['hitside_col_max'] = 0
        if track['row'].iloc[0].min() == 0:
            track['hitside_row_min'] = 1
        else:
            track['hitside_row_min'] = 0
        if track['row'].iloc[0].max() == 335:
            track['hitside_row_max'] = 1
        else:
            track['hitside_row_max'] = 0
        if track['BCID'].iloc[0].max() >= 99:
            track['hitside_BCID_max'] = 1
        else:
            track['hitside_BCID_max'] = 0

        neutron_cut = "track_energy < (0.5*length-75) and track_energy > (0.04*length-65) and track_energy > 100 and hitside_col_min == 0 and hitside_col_max == 0 and hitside_row_min == 0 and hitside_row_max == 0"

        alpha_cut = "hitside_col_min == 1 and hitside_col_max == 1 and hitside_row_min == 0 and hitside_row_max == 0 and track_energy > 800"

        xray_cut = "hitside_col_min == 0 and hitside_col_max == 0 and hitside_row_min == 0 and hitside_row_max == 0 and track_energy < 50 and length > 3000"

        if len(track.query(neutron_cut)) == 1:
            track['is_tight_neutron'] = 1
        else:
            track['is_tight_neutron'] = 0
        if len(track.query(alpha_cut)) == 1:
            track['is_tight_alpha'] = 1
        else:
            track['is_tight_alpha'] = 0
        if len(track.query(xray_cut)) == 1:
            track['is_tight_xray'] = 1
        else:
            track['is_tight_xray'] = 0
        if (len(track.query(neutron_cut)) == 0) and (len(track.query(alpha_cut)) == 0) and (len(track.query(xray_cut)) == 0):
            track['is_not_identified'] = 1
        else:
            track['is_not_identified'] = 0

        return track
        
