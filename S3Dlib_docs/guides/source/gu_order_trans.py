import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. OOPs, Order of Transforms

# 1. Define function to examine ....................................

def getCube() :  
    v = [ 
        [  0, 0, 0 ],  [  0, 1, 0 ],  [ 1 , 1, 0 ],  [  1, 0, 0 ],
        [  0, 0, 1 ],  [  0, 1, 1 ],  [ 1 , 1, 1 ],  [  1, 0, 1 ]  ]
    f = [ [0,1,2,3], [3,2,6,7], [2,1,5,6], [1,0,4,5], [0,3,7,4], [4,7,6,5] ]
    vertexCoor = np.array(v).astype(float)
    faceIndices = np.array(f)
    facecolors = np.array( ['b', 'm', 'c', 'g', 'r', 'y' ] )
    surface = s3d.Surface3DCollection(vertexCoor, faceIndices)
    surface.set_facecolor(facecolors)
    surface.transform(scale=2,translate=[-1,-1,-1])
    surface.set_surface_alpha(.3)
    return surface

scale, rot = [1,1.5,.5], s3d.eulerRot(55,20,10)

# Figure 1: ========================================================

# 2. Setup and map surface .........................................
box = [None]*4
for i in range(4) : box[i] = getCube()

box[0].transform(scale=scale)
box[1].transform(rotate=rot)
box[2].transform(scale=scale).transform(rotate=rot)
box[3].transform(rotate=rot).transform(scale=scale)

# 3. Construct figure, add surface plot ............................
title = [None]*4
title[0] = 'scale' + '\nRectangular cuboid'
title[1] = 'rotate' + '\nCube'
title[2] = r'scale$\Rightarrow$rotate'+ '\nRectangular cuboid'
title[3] = r'rotate$\Rightarrow$scale' + '\nParallelepiped'

fig = plt.figure(figsize=(4.25,4.8))
minmax = (-1.5,1.5)
for i in range(4) :
    ax = fig.add_subplot(2,2,i+1,projection='3d')
    s3d.standardAxis( ax, length=2.0 )
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
    ax.set_title(title[i])
    ax.set_proj_type('ortho')   
    ax.add_collection(box[i])
    if i == 1 : ax.add_collection3d(box[i].get_transformAxis(2,1))
    if i == 2 : ax.add_collection3d(box[i].get_transformAxis(2,1))
fig.tight_layout()

# Figure 2: ========================================================

# 2. Setup and map surface .........................................
crate = [None]*2

crate[0] = getCube().transform(rotate=rot, scale=scale)
crate[1] = getCube().transform(rotate=rot*scale)

# 3. Construct figure, add surface plot ............................
title = [None]*4
title[0] = 'transform( rotate, scale )' + '\nRectangular cuboid'
title[1] = 'transform( rotate*scale )' + '\nParallelepiped'

fig = plt.figure(figsize=(4.25,2.5))
minmax = (-1.5,1.5)
for i in range(2) :
    ax = fig.add_subplot(1,2,i+1,projection='3d')
    s3d.standardAxis( ax, length=2.0 )
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
    ax.set_title(title[i], fontsize='small')
    ax.set_proj_type('ortho')   
    ax.add_collection(crate[i])
    ax.add_collection3d(crate[i].get_transformAxis(2,1))
fig.tight_layout()

# ==================================================================
plt.show()