
import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Matplotlib Examples: Datagrid Resolution

# 1. Define function to examine .....................................

z=np.load('data/jacksboro_fault_dem.npz')['elevation']
datagrid = np.flip( z[5:50, 5:50], 0 )

# 2 & 3. Setup surfaces and plot ....................................
cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboard')
ls = s3d.elev_azim_2vector(90,-135)

fig = plt.figure(figsize=plt.figaspect(0.5))

for i in range(1,7) :
    ax = fig.add_subplot(2,3,i, projection='3d')
    ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(0,1) )
    surface = s3d.PlanarSurface(i,cmap='cardboard')
    surface.map_geom_from_datagrid( datagrid )
    surface.map_cmap_from_normals(direction=ls)
    nfaces = ' (faces: '+str(len(surface.fvIndices))+')'
    ax.set_title('rez: '+str(i)+nfaces)
    ax.set_xticks([-1,0,1])
    ax.set_yticks([-1,0,1])
    ax.set_zticks([0,1])
    ax.tick_params(labelcolor='w')

    ax.add_collection3d(surface)
    
fig.tight_layout()
plt.show()