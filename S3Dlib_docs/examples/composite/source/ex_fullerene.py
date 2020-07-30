import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Face Center Translation

# Setup surface ................................................
rez=2
size = 0.17

posSurface = s3d.SphericalSurface(basetype='dodeca')
atomPos = np.transpose(posSurface.facecenters)
interior = s3d.SphericalSurface(3,color=[0,0,0,0.5], linewidth=0).transform(scale=(1-1.6*size))
ball = None
for pos in atomPos :
    atom = s3d.SphericalSurface(rez, facecolor='lightsteelblue')
    atom.transform(scale=size, translate=pos).shade(direction=[1,1,1])
    if ball is None: ball =atom
    else:            ball += atom
total = interior+ball
total.set_edgecolor([0,0,0,0])

# Construct figure, add surface, plot ..........................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
minmax= (-.8,0.8)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()
fig.text(0.975,0.975,str(posSurface), ha='right', va='top', fontsize='smaller', multialignment='right')

ax.add_collection3d(total)

ax.view_init(elev=35, azim=-75)
plt.tight_layout()
plt.show()