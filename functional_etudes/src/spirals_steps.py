import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

# Sequence of functional steps to create the spirals figure.
#===========================================================

# 1. Define functions to examine ....................................
elev, azim = 50,-25
light_direction = [0,-1,1 ]
illum = s3d.rtv(light_direction,elev,azim)
negillum = -1.0*np.array(illum)
scolor = 'lightcoral'

def twisted_torus(rtz) :
    twists=0
    r,t,z = rtz
    ratio = .7
    phi =t*twists
    Z = ratio*np.sin(z*np.pi+phi)
    R = r + ratio*np.cos(z*np.pi+phi)
    return R,t,Z

def twisted_torus_A(rtz) :
    twists, width, radius = 4, 0.125, 0.7
    r,t,z = rtz
    ratio = radius # <--- variable assignment
    phi =t*twists
    z = width*z  # <----- added
    Z = ratio*np.sin(z*np.pi+phi)
    R = r + ratio*np.cos(z*np.pi+phi)
    return R,t,Z

def twisted_torus_B(rtz) :
    twists, width, radMax, radMin = 4, 0.125, 0.7, 0.05
    r,t,z = rtz
    #ratio = radius <---- redefined as f(t)
    ratio = radMax - t*(radMax-radMin)/(2.0*np.pi)
    phi =t*twists
    z = width*z
    Z = ratio*np.sin(z*np.pi+phi)
    R = r + ratio*np.cos(z*np.pi+phi)
    return R,t,Z

def twisted_torus_C(rtz) :
    twists, width, radMax, radMin, stretch = 4, 0.125, 0.7, 0.05, 1.4
    r,t,z = rtz
    ratio = radMax - t*(radMax-radMin)/(2.0*np.pi)
    phi =t*twists
    z = width*z
    Z = ratio*np.sin(z*np.pi+phi)
    R = r + ratio*np.cos(z*np.pi+phi)
    T = t*stretch # <--- now, extend domain
    return R,T,Z

def twisted_torus_D(rtz, rotate) :
    twists, width, radMax, radMin, stretch = 4, 0.125, 0.7, 0.05, 1.4
    r,t,z = rtz
    ratio = radMax - t*(radMax-radMin)/(2.0*np.pi)
    phi =t*twists + rotate*2*np.pi + np.pi/4.0   # <--- now f(rotate)
    z = width*z
    Z = ratio*np.sin(z*np.pi+phi)
    R = r + ratio*np.cos(z*np.pi+phi)
    T = t*stretch
    return R,T,Z

def twisted_torus_E(rtz, rotate, isInner=False) :
    twists, width, rMax, rdMin, stretch = 4, 0.125, 0.7, 0.05, 1.4
    r,t,z = rtz
    radMax, radMin = rMax, rdMin
    if isInner : 
        radMax -= 0.01
        radMin -= 0.01
    ratio = radMax - t*(radMax-radMin)/(2.0*np.pi)
    phi =t*twists + rotate*2*np.pi + np.pi/4.0
    z = width*z
    Z = ratio*np.sin(z*np.pi+phi)
    R = r + ratio*np.cos(z*np.pi+phi)
    T = t*stretch
    return R,T,Z

def Spiralset(rez, inner=False) :
    btype = 'squ_s'
    t1 = s3d.CylindricalSurface(rez, basetype=btype, color=scolor )
    t2 = s3d.CylindricalSurface(rez, basetype=btype, color=scolor )
    t3 = s3d.CylindricalSurface(rez, basetype=btype, color=scolor )
    t4 = s3d.CylindricalSurface(rez, basetype=btype, color=scolor )
    t1.map_geom_from_op( lambda rtz : twisted_torus_E(rtz,0.00,inner) )
    t2.map_geom_from_op( lambda rtz : twisted_torus_E(rtz,0.25,inner) )
    t3.map_geom_from_op( lambda rtz : twisted_torus_E(rtz,0.50,inner) )
    t4.map_geom_from_op( lambda rtz : twisted_torus_E(rtz,0.75,inner) )
    surface = t1 + t2 + t3 + t4
    return surface

# 2. Setup and map surfaces .........................................
rez=3

spiral_A = s3d.CylindricalSurface(rez, basetype='squ_s', color=scolor)
spiral_A.map_geom_from_op( lambda rtz : twisted_torus_A(rtz) )

spiral_B = s3d.CylindricalSurface(rez, basetype='squ_s', color=scolor)
spiral_B.map_geom_from_op( lambda rtz : twisted_torus_B(rtz) )

spiral_C = s3d.CylindricalSurface(rez, basetype='squ_s', color=scolor)
spiral_C.map_geom_from_op( lambda rtz : twisted_torus_C(rtz) )

spiral_D = s3d.CylindricalSurface(rez, basetype='squ_s', color=scolor)
spiral_D.map_geom_from_op( lambda rtz : twisted_torus_D(rtz,0) )

spiral_T = s3d.CylindricalSurface(rez, basetype='squ_s', color=scolor)
spiral_T.map_geom_from_op( lambda rtz : twisted_torus_D(rtz,0) )
spiral_T.shade(direction=illum).hilite(.5,direction=illum)
spiral_E = s3d.CylindricalSurface(rez, basetype='squ_s', color=scolor)
spiral_E.map_geom_from_op( lambda rtz : twisted_torus_E(rtz,0,True) )
spiral_E.shade(direction=negillum).hilite(.5,direction=negillum)
spiral_E += spiral_T

surface1 = Spiralset(rez)
surface1.shade(direction=illum).hilite(.5,direction=illum)
surface2 = Spiralset(rez,True)
surface2.shade(direction=negillum).hilite(.5,direction=negillum)
spiral_F = surface1 + surface2

grid = s3d.CylindricalSurface(5, basetype='squ_s', color=scolor)
grid.map_geom_from_op(twisted_torus)
grid.shade(direction=illum).hilite(.5,direction=illum)
grid.set_surface_alpha(.3)

surfArr = [ spiral_A, spiral_B, spiral_C, spiral_D, spiral_E, spiral_F ]
titstrg = [ 'thin torus twist',       r'radius reduction with $\theta$',
            r'extend $\theta$ domain', 'torus centered rotation',
            'add inner surface',       'duplicate with rotations']

# 3. Construct figure, add surfaces, and plot ......................
fig =plt.figure(figsize=(7,4))
maxmin = (-1,1)
for i in range(0,len(surfArr)) :
    surface = surfArr[i] + grid
    surface.set_linewidth(0)
    if i < len(surfArr)-2 :
        surface.shade(direction=illum).hilite(.5,direction=illum)
    ax = fig.add_subplot( 231+i,projection='3d')
    ax.set(xlim=maxmin, ylim=maxmin, zlim=maxmin )
    ax.view_init(elev, azim)
    ax.set_title(str(i+1)+'. '+ titstrg[i], fontsize='x-small')
    ax.set_xticks([-1,0,1])
    ax.set_yticks([-1,0,1])
    ax.set_zticks([-1,0,1])
    ax.tick_params(labelcolor='w')
    ax.add_collection3d(surface)

fig.tight_layout()
plt.show()

