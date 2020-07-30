
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from mpl_toolkits.mplot3d import axes3d
import s3dlib.surface as s3d

#.. Datagrid Order Comparison Plot

# 1. Define function to examine .....................................

X, Y, Z = axes3d.get_test_data()
Z_arr =  [ Z, np.flip(Z,0),  np.flip(Z,1), np.transpose(Z) ]
titles = [ 'default, Z', 'np.flip(Z,0)', 'np.flip(Z,1)', 'np.transpose(Z)' ]

# 2 & 3. Construct surfaces, figure, add surface, plot ..............
rez=6

fig = plt.figure(figsize=(7,6))
for i in range(4) :
    surface = s3d.PlanarSurface(rez, color='peru', basetype='oct1')
    surface.map_geom_from_datagrid( Z_arr[i] ).shade().hilite(.5)
    base = s3d.PlanarSurface(rez, cmap='Spectral', basetype='oct1')
    base.map_cmap_from_datagrid( Z_arr[i] )
    ax = fig.add_subplot(2,2,i+1,projection='3d')
    ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(0,1) )
    ax.xaxis.set_major_locator(LinearLocator(3))
    ax.yaxis.set_major_locator(LinearLocator(3))
    ax.zaxis.set_major_locator(LinearLocator(2))
    ax.set_proj_type('ortho')
    ax.set_title( titles[i], fontsize='large')
    ax.add_collection3d(surface)
    ax.add_collection3d(base)

fig.tight_layout()
plt.show()