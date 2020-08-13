import numpy as np
import pandas as pd
import root_pandas as rp
import matplotlib.pyplot as plt
import matplotlib as mpl

class event_viewer():
    def __init__(self):
        pass
    def read_file(self, inputdir = '~/data/phase3/spring_2020/maintenance_day_test/', tpc = 'tako'):
        f = inputdir + tpc + "_all.root"
        df = rp.read_root(f)
        return df
    def get_event(self, eventtype = 'neutron', entry = 0):
        data = self.read_file()
        if eventtype == 'neutron':
            data = data.loc[(data['track_energy']>50) & (data['track_energy']<60) & (data['hitside_row_min'] == 0) & (data['hitside_row_max'] == 0) & (data['hitside_col_min'] == 0) & (data['hitside_col_max'] == 0)][['column', 'row', 'BCID', 'tot']].iloc[entry]
        elif eventtype == 'alpha':
            data = data.loc[(data['hitside_row_min'] == 0) & (data['hitside_row_max'] == 0) & (data['hitside_col_min'] == 1) & (data['hitside_col_max'] == 1) & (data['theta'] <95) & (data['theta'] >85) & (data['phi'] <5) & (data['phi'] >-5)][['column', 'row', 'BCID', 'tot']].iloc[entry]
        elif eventtype == 'xray':
            data = data.loc[(data['length']>1500) & (data['track_energy']<5) & (data['hitside_row_min'] == 0) & (data['hitside_row_max'] == 0) & (data['hitside_col_min'] == 0) & (data['hitside_col_max'] == 0)][['column', 'row', 'BCID', 'tot']].iloc[entry]
        else:
            print("must enter 'neutron', 'alpha', or 'xray'. Alphas take a long time to run")
            raise ValueError
        return data
    
    def make_cubes(self, eventtype='neutron', entry = 0):
        data = self.get_event(eventtype, entry)
        xrange = data['column'].max()-data['column'].min() + 1
        yrange = data['row'].max()-data['row'].min() + 1
        zrange = data['BCID'].max()-data['BCID'].min() + 1
        x, y, z = np.indices((int(xrange), int(yrange), int(zrange)))
        points = np.concatenate((data.T['column'][:, np.newaxis], data.T['row'][:, np.newaxis], data.T['BCID'][:, np.newaxis]), axis=1)
        cubes = []
        n = len(points)
        for p in points:
            xmin = p[0]-data['column'].min()-.5
            xmax = p[0]-data['column'].min()+.5
            ymin = p[1]-data['row'].min()-.5
            ymax = p[1]-data['row'].min()+.5
            zmin = p[2]-data['BCID'].min()-.5
            zmax = p[2]-data['BCID'].min()+.5
            cubes.append((x <= xmax) & (y <= ymax) & (z <= zmax) & (x >= xmin) & (y >= ymin) & (z >= zmin))
        return cubes

    def plot(self, eventtype='neutron', entry = 0):
        data = self.get_event(eventtype, entry)
        cubes = self.make_cubes(eventtype, entry)
        xrange = data['column'].max()-data['column'].min() + 1
        yrange = data['row'].max()-data['row'].min() + 1
        zrange = data['BCID'].max()-data['BCID'].min() + 1
        fig = plt.figure(figsize=(12,8))
        ax = fig.gca(projection='3d')
        cmap = plt.cm.jet  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
        bounds = np.linspace(0, 13, 14)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        colors = []
        for val in data['tot']:
            colors.append(cmaplist[val*18])
        for i in range(0,len(cubes)):
            ax.voxels(cubes[i], facecolors = colors[i], edgecolor = 'k', linewidth = 0.3)
            print("Generated voxel %s out of %s"%(i, len(cubes)))
        ax.set_xlim(-10, xrange+10)
        ax.set_ylim(-10, yrange+10)
        ax.set_zlim(-2, zrange+2)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.set_xlabel('col')
        ax.set_ylabel('row')
        ax.set_zlabel('BCID')
        p = ax.scatter3D(data['column'], data['row'], data['BCID'], c=np.linspace(0,13,len(data['row'])), cmap = 'jet', s=0, label = 'ToT') #invisible point to map colorbar to
        plt.colorbar(p, cmap = 'jet', norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i').set_label('ToT', rotation = 270, labelpad = 20)
        plt.show()

    def plot_all_together(self): #test method using specific entries
        alpha = self.get_event("alpha", 2)
        neutron = self.get_event("neutron",19)
        xray = self.get_event("xray", 17)
        data = pd.DataFrame()
        for col in ['column', 'row', 'BCID', 'tot']:
            data[col] = list(xray[col]) + list(neutron[col]) + list(alpha[col])
        xrange = data['column'].max()-data['column'].min() + 1
        yrange = data['row'].max()-data['row'].min() + 1
        zrange = data['BCID'].max()-data['BCID'].min() + 1
        x, y, z = np.indices((int(xrange), int(yrange), int(zrange)))
        points = np.concatenate((data['column'][:, np.newaxis], data['row'][:, np.newaxis], data['BCID'][:, np.newaxis]), axis=1)
        cubes = []
        n = len(points)
        for p in points:
            xmin = p[0]-data['column'].min()-.5
            xmax = p[0]-data['column'].min()+.5
            ymin = p[1]-data['row'].min()-.5
            ymax = p[1]-data['row'].min()+.5
            zmin = p[2]-data['BCID'].min()-.5
            zmax = p[2]-data['BCID'].min()+.5
            cubes.append((x <= xmax) & (y <= ymax) & (z <= zmax) & (x >= xmin) & (y >= ymin) & (z >= zmin))
        fig = plt.figure(figsize=(12,8))
        ax = fig.gca(projection='3d')
        cmap = plt.cm.jet  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
        bounds = np.linspace(0, 13, 14)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        colors = []
        for val in data['tot']:
            colors.append(cmaplist[val*18])
        for i in range(0,len(cubes)):
            ax.voxels(cubes[i], facecolors = colors[i], edgecolor = 'k', linewidth = 0.3)
            print("Generated voxel %s out of %s"%(i, len(cubes)))
        ax.set_xlim(0, xrange)
        ax.set_ylim(-5, yrange+5)
        ax.set_zlim(-5, zrange+5)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.set_xlabel('col')
        ax.set_ylabel('row')
        ax.set_zlabel('BCID')
        p = ax.scatter3D(data['column'], data['row'], data['BCID'], c=np.linspace(0,13,len(data['row'])), cmap = 'jet', s=0, label = 'ToT') #invisible point to map colorbar to
        plt.colorbar(p, cmap = 'jet', norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i').set_label('ToT', rotation = 270, labelpad = 20)
        plt.show()
            
viewer = event_viewer()
#viewer.plot("neutron", 0)
viewer.plot_all_together()         
