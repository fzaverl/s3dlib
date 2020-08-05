import copy
import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Guides: Shading, Highlighting and Color Mapped Normals -
#   Color Mapping Normals: Default, Cmaps, Shade v. Cmap

rez = 4
minmax = (-.7,.7)
height = 1.6

#. Figure 1 - default ...............................................
rez = 4
surface = s3d.SphericalSurface(rez)
surface.map_cmap_from_normals(direction= [1,0.2,1] )

fig = plt.figure(figsize=(height ,height))
ax = plt.axes( projection='3d')
ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax.set_title('Default (viridis)')
ax.set_axis_off()

ax.add_collection3d(surface)

fig.tight_layout()

#. Figure 2 - cmaps .................................................
cmaps = ['plasma', 'magma', 'winter', 'autumn', 'hsv' ]


fig = plt.figure(figsize=(6.75,height))
for i in range( len(cmaps) ) :
    surface = s3d.SphericalSurface(rez,cmap=cmaps[i])
    surface.map_cmap_from_normals(direction=[1,0.2,1])
    ax = fig.add_subplot( 151+i , projection='3d')
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    ax.set_title(str(cmaps[i]) , fontsize='large')
    ax.set_axis_off()

    ax.add_collection3d(surface)

fig.tight_layout()

#. Figure 3 - Shade v. Cmap .........................................
rez = 4
dark, light  = [0.25,0.15,0], [1,.9,.75]
cmu.rgb_cmap_gradient(dark,light,'cardboard')

shaded_surface = s3d.SphericalSurface(rez,color=light)
shaded_surface.shade(0.25,direction= [1,0.2,1] )
cmapped_surface = s3d.SphericalSurface(rez,cmap='cardboard')
cmapped_surface.map_cmap_from_normals(direction= [1,0.2,1] )

fig = plt.figure(figsize=(2.8 ,height))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')
ax1.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax2.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax1.set_title('Shade')
ax2.set_title('Cmap')
ax1.set_axis_off()
ax2.set_axis_off()

ax1.add_collection3d(shaded_surface)
ax2.add_collection3d(cmapped_surface)

fig.tight_layout()

#....................................................................
plt.show()