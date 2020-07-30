import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. High Resolution Surface Plot

# 1. Define function to examine .....................................

def mayDemo(rtp) :
    r,theta,phi = rtp
    m0 = 4; m1 = 3; m2 = 2; m3 = 3; m4 = 6; m5 = 2; m6 = 6; m7 = 4;
    R = np.sin(m0*phi)**m1 + np.cos(m2*phi)**m3 + np.sin(m4*theta)**m5 + np.cos(m6*theta)**m7
    return R,theta,phi

def yaxisDir(rtp) :
    x,y,z = s3d.SphericalSurface.coor_convert(rtp,True)
    return y

# 2. Setup and map surfaces .........................................
rez = 8
cmap = cmu.hue_cmap(2.0,'r','b')

surface = s3d.SphericalSurface(rez, basetype='dodeca')
surface.map_geom_from_op(mayDemo)
surface.map_cmap_from_op(yaxisDir,cmap).shade(.2).hilite(.7,focus=2)
surface.transform(rotate=s3d.eulerRot(0,-90,False))

# 3. Construct figures, add surface, plot ...........................

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975,str(surface), ha='right', va='top',
        fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
ax.set(xlim=(-1.2,1.2), ylim=(-0.8,1.6), zlim=(-1.2,1.2))
ax.view_init(30,45)
ax.set_axis_off()

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()