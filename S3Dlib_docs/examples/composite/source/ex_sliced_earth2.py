import copy
import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Different Sub-surface Type, clip & 5 surfaces.

# 1. Define functions to examine ....................................

def addFace(rtp) :
    r,t,p = rtp
    toAdd = np.logical_and(t>3*np.pi/2, p<np.pi/2 )
    return toAdd

def revFace(rtp) :
    rev = np.logical_not(addFace(rtp))
    return rev

# 2. Setup and map surfaces .........................................
rez=7

exterior = s3d.SphericalSurface(rez, basetype='octa')
exterior.map_color_from_image('data/blue_marble.png')
exterior.clip( revFace ).shade(.5,[0,0,1])

interior1 = s3d.PolarSurface(5, basetype='squ', cmap="hot")
interior1.map_cmap_from_op( lambda rtz : 1-rtz[0]  )
interior2 = copy.copy(interior1).transform( [ [0,0,-1], [0,1,0], [1,0,0] ] )
interior3 = copy.copy(interior1).transform( [ [1,0, 0], [0,0,1], [0,1,0] ] )
interior = interior1 + interior2 + interior3

cmap = cmu.hsv_cmap_gradient('orange','gold')
core = s3d.SphericalSurface(rez, basetype='octa', color='gold').transform(scale=.325)
core.map_cmap_from_normals(cmap).hilite(.5,[1,-1,1])

surface = interior + exterior + core

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1), facecolor='black' )
desc = str(surface) + '\n' + str(exterior) + '\n' + str(interior)
fig.text(0.975,0.975, desc, ha='right', va='top', 
        fontsize='smaller', multialignment='right', color='white')
ax = fig.add_subplot(111, projection='3d')
minmax = (-0.75,0.75)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()
ax.set_facecolor('black')

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()