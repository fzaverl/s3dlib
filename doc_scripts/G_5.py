import copy
import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Guides: Shading, Highlighting and Color Mapped Normals -
#   Highlighting: Height and Focus

rez = 4
minmax = (-.7,.7)

#. Figure 1 - Height ................................................
heightLevel = [0, 0.25, 0.5, 0.75, 1 ]

surface_ref = s3d.SphericalSurface(rez,color='peru')
surface_ref.map_color_from_image('data/earth.png')

width,height = 6.75, 1.6
fig = plt.figure(figsize=(width,height))
for i in range( len(heightLevel) ) :
    surface = copy.copy(surface_ref)
    surface.hilite( heightLevel[i], direction=[1,0.2,1])
    ax = fig.add_subplot( 151+i , projection='3d')
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    ax.set_title(str(heightLevel[i]) , fontsize='large')
    ax.set_axis_off()

    ax.add_collection3d(surface)

fig.tight_layout()

#. Figure 2 - Focus .................................................
focusLevel = [1,1.25,1.67,2,  1, 0.67, 0.5, 0.25 ]

width,height = 5.25, 3.0
fig = plt.figure(figsize=(width,height))
for i in range( len(focusLevel) ) :
    surface = s3d.SphericalSurface(rez,color='peru')
    surface.hilite( direction=[1,0.2,1], focus=focusLevel[i])
    ax = fig.add_subplot( 241+i , projection='3d')
    ax.set_title(str(focusLevel[i]) , fontsize='large')
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    ax.set_axis_off()

    ax.add_collection3d(surface)

fig.tight_layout()

#....................................................................
plt.show()