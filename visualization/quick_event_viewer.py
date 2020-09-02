# based of of script https://stackoverflow.com/questions/42611342/representing-voxels-with-matplotlib
import root_pandas as rp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

def plot_event(data, track_type = 'alpha', evt_num = 0, cmap = plt.cm.jet):
    #define generic cuts
    neutron = data.loc[(data['hitside_row_min'] == 0) & (data['hitside_row_max'] == 0) & (data['hitside_col_min'] == 0) & (data['hitside_col_max'] == 0) & (data['track_energy'] > 50)][['column', 'row', 'BCID','tot']] #selects a nice looking nuclear recoil

    xray = data.loc[(data['hitside_row_min'] == 0) & (data['hitside_row_max'] == 0) & (data['hitside_col_min'] == 0) & (data['hitside_col_max'] == 0) & (data['track_energy'] < 50) & (data['length'] > 10000)][['column', 'row', 'BCID','tot']] #selects a long xray with large gaps

    alpha = data.loc[(data['hitside_row_min'] == 0) & (data['hitside_row_max'] == 0) & (data['hitside_col_min'] == 1) & (data['hitside_col_max'] == 1) & (data['track_energy'] > 300) & (data['theta'] > 80) & (data['theta'] < 100)][['column', 'row', 'BCID','tot']] #selects a relatively horizontal alpha

    ################################################################

    # Define track to be plotted based off of input track_type
    if track_type.lower() == 'alpha':
        plotter = alpha.iloc[evt_num]
    elif track_type.lower() == 'neutron':
        plotter = neutron.iloc[evt_num]
    else:
        plotter = xray.iloc[evt_num]

    # Define and discretize colormap to plot
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap_discrete = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = np.linspace(0, 13, 14)
    norm = mpl.colors.BoundaryNorm(bounds, cmap_discrete.N)
    colors = []
    for val in plotter[3]:
        colors.append(cmaplist[val*18])

    def cuboid_data(o, size=(1,1,1)):
        X = [[[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
         [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
         [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
         [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
         [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
         [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]]]
        X = np.array(X).astype(float)
        for i in range(3):
            X[:,:,i] *= size[i]
        X += np.array(o)
        return X

    def plotCubeAt(positions,sizes=None,colors=None,edgecolors=None, **kwargs):
        if not isinstance(colors,(list,np.ndarray)): colors=["C0"]*len(positions)
        if not isinstance(sizes,(list,np.ndarray)): sizes=[(1,1,1)]*len(positions)
        if not isinstance(edgecolors,(list,np.ndarray)): edgecolors=[[0,0,0,.5]]*len(positions)
        g = []
        for p,s,c in zip(positions,sizes,colors):
            g.append( cuboid_data(p, size=s) )
        return Poly3DCollection(np.concatenate(g), facecolors=np.repeat(colors,6, axis=0),
                                edgecolors=np.repeat(edgecolors,6, axis=0), lw=0.25, **kwargs)

    ax = plt.gca(projection='3d')
    ax.w_zaxis.line.set_lw(0.)
    ax.set_zticks([])
    ### things to plot
    x = plotter[0] * 250.
    y = ( 335 - plotter[1] ) * 50.
    z = plotter[2] * 250.
    ###
    positions = np.array([ x , y , z ]).transpose()
    ###
    pc = plotCubeAt(positions, sizes=[(250, 50, 250)]*len(positions), colors=colors)
    ax.add_collection3d(pc)
    ################################################################
    # set bounds to max axis & modify aspect ratio
    ax.set_xlim(0,20000)
    ax.set_ylim(0,16800)
    ax.set_zlim(0,25000)
        ################################################################
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    #ax.set_zticklabels([])
    return ax, x, y, z

def make_colorbar(cmap):
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap_discrete = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = np.linspace(0, 13, 14)
    norm = mpl.colors.BoundaryNorm(bounds, cmap_discrete.N)
    ax = plt.gca(projection='3d')
    col = ax.scatter3D(bounds, bounds, bounds, c=bounds, cmap = cmap, s=0, label = 'ToT')
    cbar = plt.colorbar(col, cmap = 'jet', norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i', orientation = 'horizontal').set_label('ToT', labelpad = 20)
    return cbar

#Read in data
tako = rp.read_root("~/data/phase3/spring_2020/05-09-20/tpc_root_files/tako_all.root")

# Create colormap
cmap = plt.cm.jet

#plot voxelized events
ax1, x1, y1, z1 = plot_event(tako, 'alpha', 2, cmap = cmap)
ax2, x2, y2, z2 = plot_event(tako, 'xray', 3, cmap = cmap)
ax3, x3, y3, z3 = plot_event(tako, 'neutron', 11, cmap = cmap) #return ax to scale axes

#determine max and min of each axis
xmin = np.min([np.min(x1), np.min(x2), np.min(x3)])
xmax = np.max([np.max(x1), np.max(x2), np.max(x3)])
ymin = np.min([np.min(y1), np.min(y2), np.min(y3)])
ymax = np.max([np.max(y1), np.max(y2), np.max(y3)])
zmin = np.min([np.min(z1), np.min(z2), np.min(z3)])
zmax = np.max([np.max(z1), np.max(z2), np.max(z3)])

#set minor tick marks at the pixel level
ax3.yaxis.set_major_locator(MultipleLocator(80))
ax3.xaxis.set_major_locator(MultipleLocator(336))
#ax3.xaxis.set_minor_locator(MultipleLocator(250))
#ax3.yaxis.set_minor_locator(MultipleLocator(50))

#Scale axes
#ax3.set_xlim(xmin,xmax)
#ax3.set_ylim(ymin,ymax)
#ax3.set_zlim(zmin,zmax)

cbar = make_colorbar(cmap)

plt.show()
