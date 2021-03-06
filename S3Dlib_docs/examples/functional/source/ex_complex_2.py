import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Complex Number Representation, Geometry and Colormap: Planar Coordinates

# 1. Define functions to examine ....................................

real = True
imaginary =  not real

def sqrt_Z(xyz, isReal, isPos=True) :
    x,y,z = xyz
    r = np.sqrt( x*x + y*y )
    t = np.arctan2(y,x)
    if isReal :  Z = np.sqrt(r)*np.cos(t/2)
    else :       Z = np.sqrt(r)*np.sin(t/2)
    if not isPos : Z = -Z 
    return x,y,Z

def square_Z(xyz, isReal) :
    x,y,z = xyz
    if isReal : Z = x*x - y*y
    else :      Z = 2*x*y
    return x,y,Z

# 2. Setup and map surfaces .........................................

zMap = cmu.hsv_cmap_gradient( [.5,1,1], [1.5,1,1] )

surface_1pos = s3d.PlanarSurface(3, basetype='oct1', cmap=zMap,linewidth=0.1)
surface_1pos.map_cmap_from_op( lambda xyz : sqrt_Z(xyz,imaginary)[2] )
surface_1pos.map_geom_from_op( lambda xyz : sqrt_Z(xyz,real) )

surface_1neg = s3d.PlanarSurface(3, basetype='oct1', cmap=zMap,linewidth=0.1)
surface_1neg.map_cmap_from_op( lambda xyz : sqrt_Z(xyz,imaginary,False)[2] )
surface_1neg.map_geom_from_op( lambda xyz : sqrt_Z(xyz,real,False) )

surface_1 = surface_1pos + surface_1neg
surface_1.set_linewidth(0.1)
surface_1.set_edgecolor([0,0,0])

surface_I = s3d.PlanarSurface(5, facecolor='red')
surface_I.map_geom_from_op( lambda xyz : square_Z(xyz,imaginary) )
surface_R = s3d.PlanarSurface(5, facecolor='blue')
surface_R.map_geom_from_op( lambda xyz : square_Z(xyz,real) )

surface_2 = surface_I + surface_R
surface_2.shade(direction=[1,1,1])

# 3. Construct figure, add surfaces, and plot .....................

fcc = r'     f: $\mathrm{\mathbb{C}}$ $\to$  $\mathrm{\mathbb{C}}$' +'\n\n'
minmax, ticks = (-1,1) ,  [-1,0,1]
fig = plt.figure(figsize=plt.figaspect(0.6/1.2))

ax1 = fig.add_subplot(121, projection='3d')
ax1.view_init(20, 205)
ax1.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax1.set_xticks(ticks)
ax1.set_yticks(ticks)
ax1.set_zticks(ticks)
ax1.set_xlabel('X axis')
ax1.set_ylabel('Y axis')
ax1.set_title( r'f(z) = $\sqrt{z}$' + fcc )
cbar = plt.colorbar(surface_1pos.cBar_ScalarMappable, ax=ax1,  shrink=0.6, pad=-.05 )
cbar.set_label('Imaginary', rotation=270, labelpad = 15)
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%5.2f'))
real_handle = mlines.Line2D([], [], color='grey', marker='+', linestyle='None',
                          markersize=7, label='Real')
ax1.legend(handles=[real_handle])

ax2 = fig.add_subplot(122, projection='3d')
ax2.set(xlim=minmax, ylim=minmax, zlim=(-2,2))
ax2.set_xticks(ticks)
ax2.set_yticks(ticks)
ax2.set_zticks([-2,-1,0,1,2])
ax2.set_xlabel('X axis')
ax2.set_ylabel('Y axis')
ax2.set_title( r'f(z) =  $z^2$' + fcc )
red_patch = mpatches.Patch(color='red', label='Img')
blue_patch = mpatches.Patch(color='blue', label='Real')
ax2.legend(handles=[red_patch,blue_patch])

ax1.add_collection3d(surface_1)
ax2.add_collection3d(surface_2)

fig.tight_layout()
plt.show()