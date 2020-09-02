# based of of script https://stackoverflow.com/questions/42611342/representing-voxels-with-matplotlib
import root_pandas as rp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

tako = rp.read_root("~/data/phase3/spring_2020/05-09-20/tpc_root_files/tako_all.root")

neutron = tako.loc[(tako['hitside_row_min'] == 0) & (tako['hitside_row_max'] == 0) & (tako['hitside_col_min'] == 0) & (tako['hitside_col_max'] == 0) & (tako['track_energy'] > 300)][['column', 'row', 'BCID','tot']].iloc[1] #selects a nice looking nuclear recoil

xray = tako.loc[(tako['length']>1500) & (tako['track_energy']<5) & (tako['hitside_row_min'] == 0) & (tako['hitside_row_max'] == 0) & (tako['hitside_col_min'] == 0) & (tako['hitside_col_max'] == 0)][['column', 'row', 'BCID','tot']].iloc[0] #selects a long xray with large gaps

alpha = tako.loc[(tako['hitside_row_min'] == 0) & (tako['hitside_row_max'] == 0) & (tako['hitside_col_min'] == 1) & (tako['hitside_col_max'] == 1) & (tako['phi'] <5) & (tako['phi'] >-5) & (tako['track_energy'] > 700)][['column', 'row', 'BCID','tot']].iloc[1] #selects a relatively horizontal alpha

tester = tako.loc[(tako['hitside_BCID_max'] == 1)][['column', 'row', 'BCID','tot']].iloc[1]
################################################################
# track being plotted
plotter = alpha
################################################################
# these are the rgb values being used
'''
rgb color gradient from blue (0) -> red (13)

0 0000FF
1 004EFF
2 009CFF
3 00EBFF
4 00FFC4
5 00FF75
6 00FF27
7 27FF00
8 75FF00
9 C4FF00
10 FFEB00
11 FF9C00
12 FF4E00
13 FF0000
'''
################################################################
# dictionaries that contain the decimal values of the color scale

cmap = plt.cm.jet
cmaplist = [cmap(i) for i in range(cmap.N)]
cmap_discrete = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
bounds = np.linspace(0, 13, 14)
norm = mpl.colors.BoundaryNorm(bounds, cmap_discrete.N)
colors = []
for val in plotter[3]:
    colors.append(cmaplist[val*18])

edges_dict = {
    0 : np.array([0,0,1,1]),
    1: np.array([0,.31,1,1]),
    2: np.array([0,.61,1,1]),
    3: np.array([0,.92,1,1]),
    4: np.array([0, 1, .77,1]),
    5: np.array([0, 1, .46,1]),
    6: np.array([0, 1, .15,1]),
    7: np.array([.15, 1, 0,1]),
    8: np.array([.46, 1, 0,1]),
    9: np.array([.77, 1, 0,1]),
    10: np.array([1, .92, 0,1]),
    11: np.array([1, .61, 0,1]),
    12: np.array([1, .31, 0,1]),
    13: np.array([1, 0, 0,1])
}
# replace tot values with rbg list
edgecolors = pd.Series(plotter[3].transpose()).map(edges_dict).to_numpy()
################################################################
# this part is mainly from stackexchange
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

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
    return Poly3DCollection(np.concatenate(g),  
                            facecolors=np.repeat(colors,6, axis=0), edgecolors=np.repeat(edgecolors,6, axis=0), lw=0.25, **kwargs)

fig = plt.figure()
ax = fig.gca(projection='3d', proj_type='ortho')
#ax.set_aspect('equal')#broken?
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
x_diff = np.max(x) - np.min(x)
x_mid = ( np.min(x) + np.max(x) ) / 2
y_diff = np.max(y) - np.min(y)
y_mid = ( np.min(y) + np.max(y) ) / 2
z_diff = np.max(z) - np.min(z)
z_mid = ( np.min(z) + np.max(z) ) / 2
maxi = np.max([x_diff, y_diff, z_diff])
if maxi == x_diff:
    ax.set_xlim([np.min(x), np.max(x)])
    ax.set_ylim([y_mid - x_diff / 2, y_mid + x_diff / 2])
    ax.set_zlim([z_mid - x_diff / 2, z_mid + x_diff / 2])
elif maxi == y_diff:
    ax.set_ylim([np.min(y), np.max(y)])
    ax.set_xlim([x_mid - y_diff / 2, x_mid + y_diff / 2])
    ax.set_zlim([z_mid - y_diff / 2, z_mid + y_diff / 2])
else:
    ax.set_zlim([np.min(z), np.max(z)])
    ax.set_xlim([x_mid - z_diff / 2, x_mid + z_diff / 2])
    ax.set_ylim([y_mid - z_diff / 2, y_mid + z_diff / 2])   
################################################################
ax.set_xticklabels([i for i in range(0,80)])
ax.set_yticklabels([i for i in range(0,336)])
ax.set_zticklabels([])
################################################################
col = ax.scatter3D(x, y, z, c=np.linspace(0,13,len(x)), cmap = cmap, s=0, label = 'ToT')
plt.colorbar(col, cmap = 'jet', norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i', orientation = 'horizontal').set_label('ToT', labelpad = 20)
plt.show()
