import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Alpha Channel Adjustments

# 1. Define functions to examine ....................................

def ringDef(rtz) :
    r,t,z = rtz
    minRad, maxRad = 1.28, 2.41  # normalized to saturn radius.
    m = (maxRad-minRad)/2.0
    b = (maxRad+minRad)/2.0
    R = m*z + b
    Z = np.zeros(len(z))
    return R,t,Z

# 2. Setup and map surfaces .........................................
rez=4

surface = s3d.SphericalSurface(rez)
surface.map_color_from_image('data/saturn_surface.png')
surface.shade(direction=[1,1,1])
surfaceInfo = str(surface)

ring = s3d.CylindricalSurface(rez+2)
ring.map_color_from_image('data/saturn_rings_trans.png')
ring.set_surface_alpha(0.1)
ring.map_geom_from_op(ringDef)
ringInfo = str(ring)

saturn = surface + ring
saturn.transform(rotate=s3d.eulerRot(0,30))
info = str(saturn) + '\n' + surfaceInfo + '\n' + ringInfo

# 3. Construct figures, add surfaces, and plot .....................

fig = plt.figure(figsize=plt.figaspect(1.0), facecolor='black' )
fig.text(0.975,0.975,info, ha='right', va='top',
        fontsize='smaller', multialignment='right', color='white')
ax = fig.add_subplot(111, projection='3d')
ax.view_init(0, -70)
minmax = (-1.5,1.5)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()
ax.set_facecolor('black')

ax.add_collection3d(saturn)

fig.tight_layout()
plt.show()