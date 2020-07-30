import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Clipping a surface

# 1. Define function to examine .....................................

def addFace(rtp) :
    r,t,p = rtp
    toAdd = np.logical_and(t>3*np.pi/2, p<np.pi/2 )
    return toAdd

def revFace(rtp) :
    rev = np.logical_not(addFace(rtp))
    return rev

# 2. Setup and map surfaces .........................................
rez = 5
cmap = cmu.hsv_cmap_gradient( [0,1,1], [0.333,1,.65], smooth=1.6 )

move, amount = [1,-1,1], 0.3
moveOut = np.multiply(move,amount)
movebck = np.multiply(move,-amount)

surface = s3d.SphericalSurface(rez, basetype='octa').clip(revFace)
surface.transform(translate=movebck)

chip = s3d.SphericalSurface(rez, basetype='octa').clip(addFace)
chip.transform(translate=moveOut)

info = str(surface)+'\n'+str(chip)
surface = surface+chip
surface.map_cmap_from_normals(cmap,movebck).shade(.5,direction=[0,0,1])
info = info+'\n'+str(surface)

# 3. Construct figures, add surface, plot ...........................

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975,info, ha='right', va='top', fontsize='smaller', multialignment='right')

ax = plt.axes(projection='3d')
minmax = (-1.0,1.0)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.view_init(elev=15, azim=-75)
ax.set_xticks([-1,0,1])
ax.set_yticks([-1,0,1])
ax.set_zticks([-1,0,1])

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()
