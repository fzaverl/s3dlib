import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

#.. Base Class Surface

# 1. Define function to examine ....................................

v = [ 
    [  0, 0, 0 ],  [  0, 1, 0 ],  [ 1 , 1, 0 ],  [  1, 0, 0 ],
    [  0, 0, 1 ],  [  0, 1, 1 ],  [ 1 , 1, 1 ],  [  1, 0, 1 ]  ]
f = [ [0,1,2,3], [3,2,6,7], [2,1,5,6], [1,0,4,5], [0,3,7,4], [4,7,6,5] ]
e = [ [3,2], [2,1], [1,0], [0,3],   [7,6], [6,5], [5,4], [4,7],   [2,6], [1,5], [0,4], [3,7]  ]
vertexCoor = np.array(v).astype(float)
faceIndices = np.array(f)
edgeIndices = np.array(e)
facecolors = np.array( ['b', 'm', 'c', 'g', 'r', 'y' ] )

# 2. Setup and map surface .........................................

surface = s3d.Surface3DCollection(vertexCoor, faceIndices, edgeIndices, facecolors=facecolors)
surface.transform(scale=2,translate=[-1,-1,-1])
surface.set_surface_alpha(.5)

# 3. Construct figure, add surface plot ............................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
minmax = (-1.5,1.5)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_title(str(surface))
ax.set_proj_type('ortho')

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()