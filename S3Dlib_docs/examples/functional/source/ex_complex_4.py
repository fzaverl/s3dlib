import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Complex Number Representation, Geometry and Colormap: Alternative Represenations

# 1. Define functions to examine ....................................

real = True
imaginary =  not real

def sqrt_Z(rtz, isReal) :
    r,t,z = rtz
    T = 2*np.mod( (t + 2*np.pi) , 2*np.pi)
    if isReal :  Z = np.sqrt(r)*np.cos(T/2)
    else :       Z = np.sqrt(r)*np.sin(T/2)
    return r,T,Z
    
def square_Z(rtz, isReal) :
    r,t,z = rtz
    if isReal : Z = r*r*np.cos(2*t)
    else :      Z = r*r*np.sin(2*t)
    return r,t,Z

# 2. Setup and map surfaces .........................................

zMap = cmu.hsv_cmap_gradient( [.5,1,1], [1.5,1,1] )

surface1_I = s3d.PolarSurface(6, basetype='squ', facecolor='red')
surface1_I.map_geom_from_op( lambda rtz : sqrt_Z(rtz,imaginary) )
surface1_R = s3d.PolarSurface(6, basetype='squ', facecolor='blue')
surface1_R.map_geom_from_op( lambda rtz : sqrt_Z(rtz,real) )
surface_1 = (surface1_I + surface1_R).shade(direction=[1,1,1])

surface_2 = s3d.PolarSurface(4, cmap=zMap,linewidth=0.1)
surface_2.map_cmap_from_op( lambda rtz : square_Z(rtz,imaginary)[2] )
surface_2.map_geom_from_op( lambda rtz : square_Z(rtz,real) )
surface_2.set_edgecolor([0,0,0])

# 3. Construct figure, add surfaces, and plot .....................

fcc = r'     f: $\mathrm{\mathbb{C}}$ $\to$  $\mathrm{\mathbb{C}}$' +'\n\n'
minmax, ticks = (-1,1) ,  [-1,0,1]
fig = plt.figure(figsize=plt.figaspect(0.6/1.2))
ax1 = fig.add_subplot(121, projection='3d')

ax1.view_init(30, 205)
ax1.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax1.set_xticks(ticks)
ax1.set_yticks(ticks)
ax1.set_zticks(ticks)
ax1.set_xlabel('X axis')
ax1.set_ylabel('Y axis')
ax1.set_title( r'f(z) = $\sqrt{z}$' + fcc )
red_patch = mpatches.Patch(color='red', label='Img')
blue_patch = mpatches.Patch(color='blue', label='Real')
ax1.legend(handles=[red_patch,blue_patch])

ax2 = fig.add_subplot(122, projection='3d')
ax2.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax2.set_xticks(ticks)
ax2.set_yticks(ticks)
ax2.set_zticks(ticks)
ax2.set_xlabel('X axis')
ax2.set_ylabel('Y axis')
ax2.set_title( r'f(z) =  $z^2$' + fcc )
cbar = plt.colorbar(surface_2.cBar_ScalarMappable, ax=ax2,  shrink=0.6, pad=-.05 )
cbar.set_label('Imaginary', rotation=270, labelpad = 15)
cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%5.2f'))
real_handle = mlines.Line2D([], [], color='grey', marker='+', linestyle='None',
                          markersize=7, label='Real')
ax2.legend(handles=[real_handle])

ax1.add_collection3d(surface_1)
ax2.add_collection3d(surface_2)
fig.tight_layout()
plt.show()