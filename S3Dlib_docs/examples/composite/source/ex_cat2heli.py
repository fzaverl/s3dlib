import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Parametric set of Surfaces 2

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
rez = 4
cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboard')

surface = s3d.CylindricalSurface(rez, basetype='squ_s', antialiased=True)
surface.map_geom_from_op( lambda rtz : catenoid_helicoid(rtz,1), returnxyz=True )
surface.map_cmap_from_op( lambda rtz : rtz[0],'cardboard')

mesh = s3d.CylindricalSurface(rez, basetype='squ_s', facecolor=[0.25,0.15,0,0.2], edgecolor=[0,0,0,0])
mesh.map_geom_from_op( lambda rtz : catenoid_helicoid(rtz,0.85), returnxyz=True )

# 3. Construct figures, add surface, plot ...........................

fig = plt.figure(figsize=plt.figaspect(1))
fig.text(0.975,0.975,str(surface), ha='right', va='top', fontsize='smaller', multialignment='right')
ax = plt.axes(projection='3d')
minmax = (-1.2,1.2)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax )
ax.view_init(25, -28)
ax.set_axis_off()

ax.add_collection3d(surface + mesh)

fig.tight_layout()
plt.show()