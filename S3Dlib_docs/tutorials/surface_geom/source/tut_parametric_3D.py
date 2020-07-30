import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import s3dlib.surface as s3d

# 1. Define function to examine ....................................

def planarfunc(xyz,N) :
    x,y,z = xyz
    r = np.sqrt( x**2 + y**2)
    Z = ( np.sin(6*r) + N )/2 - 1
    return x,y,Z

# 2. Setup and map surface .........................................

y1 = s3d.PlanarSurface(4, color='C0' )
y1.map_geom_from_op( lambda xyz : planarfunc(xyz, 1) ).shade()
y2 = s3d.PlanarSurface(4, color='C1' )
y2.map_geom_from_op( lambda xyz : planarfunc(xyz, 2) ).shade()
y3 = s3d.PlanarSurface(4, color='C2')
y3.map_geom_from_op( lambda xyz : planarfunc(xyz, 3) ).shade()

# 3. Construct figure, add surface, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
minmax, ticks = (-1.2,1.2), [-1,0,1]
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_xticks( ticks )
ax.set_yticks( ticks )
ax.set_zticks( ticks )

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z = f(X,Y,N)')
ax.set_title( r'Z =  $\frac{1}{2} [ sin(6r) +  N ] - 1 $'+
              r'   $\mathit{where}\  r = \left( X^2 + Y^2 \right)^\frac{1}{2} $'+'\n\n' )

C0_patch = mpatches.Patch(color='C0', label='N = 1')
C1_patch = mpatches.Patch(color='C1', label='N = 2')
C2_patch = mpatches.Patch(color='C2', label='N = 3')
ax.legend(handles=[C0_patch,C1_patch,C2_patch])

ax.add_collection3d( y1 + y2 + y3 )

plt.show()
