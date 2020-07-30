import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Surface Displacement Vector Field

# 1. Define function to examine .....................................

def catenoid_helicoid(rtz, A) :
    r,t,z = rtz
    A = A*np.pi  #  -1 < A < 1
    cosA, sinA = np.cos(A), np.sin(A)
    U, V = t, z   
    x =  cosA * np.sinh(V) * np.sin(U) +   sinA * np.cosh(V) * np.cos(U)
    y = -cosA * np.sinh(V) * np.cos(U) +   sinA * np.cosh(V) * np.sin(U)
    Z = ( U/np.pi- 1.0 ) *cosA +  V * sinA
    return x,y,Z

# 2. Setup and map surfaces .........................................
rez = 5
cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboard')
card_trans = cmu.alpha_cmap( 'cardboard', 0.05 )

surface = s3d.CylindricalSurface(rez, basetype='squ_s')
surface.map_geom_from_op( lambda rtz : catenoid_helicoid(rtz,.5), returnxyz=True )
surface.map_cmap_from_op( lambda rtz : rtz[0], card_trans)

low_rez = s3d.CylindricalSurface(2, basetype='squ_s')
low_rez.map_geom_from_op( lambda rtz : catenoid_helicoid(rtz,.5), returnxyz=True )
vf = low_rez.dispfield_from_op(lambda rtz : catenoid_helicoid(rtz,0.6), returnxyz=True, scale=1 )

# 3. Construct figures, add surface, plot ...........................

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975,str(vf), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
minmax = (-1.2,1.2)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax.view_init(25, -28)

ax.add_collection3d(surface)
ax.add_collection3d(vf)

fig.tight_layout()
plt.show()