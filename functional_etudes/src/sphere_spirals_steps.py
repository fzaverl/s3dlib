import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d
import copy

# Sequence of functional steps to create the sphere spirals figure.
#==================================================================

# 1. Define function to examine .....................................
ttlcolor =  [0.502, 0.200, 0.278, 0.1] 
fgbgcolor = [0.957, 0.925, 0.882]
outercolor = [0.875, 0.682, 0.388]
innercolor = [0.847, 0.525, 0.471]

elev, azim = 45,30
light_direction = [1,-1,1 ]
illum = s3d.rtv(light_direction,elev,azim)
negillum = -1.0*np.array(illum)

def cliptwo(rtp) :
    r,t,p = rtp
    val = np.full(len(r),False)
    val1 = np.logical_and(t>0.0,t<np.pi/2)
    val2 = np.logical_and(t>np.pi,t<3*np.pi/2)
    val = np.logical_or(val1,val2)
    return val

def opensphere(rtp, rtt = 0.0) :
    r,t,p = rtp
    T = t/2.0 + rtt   
    return r,T,p

def twist(rtp, twist) :
    r,t,p = rtp
    f = np.cos(t) + 1
    b =  2*p/np.pi -1
    g = ( 1 + np.sign(b)*np.power(  np.abs(b), 6 ) )/2
    #T = t - 2*twist*p  .... twist rate is not llinear, so...
    T = t - 2*twists*(p+g*2*np.pi)
    return r,T,p

# 2. Setup and map surfaces .........................................
rez = 5
twists = 1

surface1 = s3d.SphericalSurface(rez, basetype='octa', color=outercolor)
surface1.clip( cliptwo )
surface_R = copy.copy(surface1).shade()

surface1.map_geom_from_op(lambda rtp : opensphere(rtp))
surface_S = copy.copy(surface1).shade()
surface1.map_geom_from_op(lambda rtp : twist(rtp,twists))
surface_T = copy.copy(surface1).shade()

surface2 = s3d.SphericalSurface(rez, basetype='octa', color=outercolor)
surface2.clip( cliptwo )
surface2.map_geom_from_op(lambda rtp : opensphere(rtp, np.pi))
surface2.map_geom_from_op(lambda rtp : twist(rtp,twists))

surfaceA = surface1 + surface2
surfaceA.shade(direction=illum).hilite(.5,direction=illum)
surface_U = copy.copy(surfaceA)

surface1 = s3d.SphericalSurface(rez, basetype='octa', color=innercolor)
surface1.clip( cliptwo )
surface1.map_geom_from_op(lambda rtp : opensphere(rtp))
surface1.map_geom_from_op(lambda rtp : twist(rtp,twists))

surface2 = s3d.SphericalSurface(rez, basetype='octa', color=innercolor)
surface2.clip( cliptwo )
surface2.map_geom_from_op(lambda rtp : opensphere(rtp, np.pi))
surface2.map_geom_from_op(lambda rtp : twist(rtp,twists))

surfaceB = surface1 + surface2
surfaceB.transform(scale=0.95)
surfaceB.shade(direction=negillum).hilite(.5,direction=negillum)

surface = surfaceA + surfaceB
surface_V = copy.copy(surface)
grid = s3d.SphericalSurface(4, basetype='icosa', color='k')
grid.map_cmap_from_normals('binary',illum)
grid.shade(direction=negillum).hilite(direction=negillum)
grid.set_surface_alpha(0.05) 
grid.transform(scale=1.01)

surface = surface + grid

surfArr = [ surface_R, surface_S, surface_T, surface_U, surface_V, surface ]
titstrg = [ 'clip two quadrants',  'map to half sphere',
            'twist',               'add other half rotated',
            'add inner surface',   'add outer grid']

# 3. Construct figure, add surface, plot ............................

fig =plt.figure(figsize=(6.5,4))
maxmin = (-.7,.7)
for i in range(0,len(surfArr)) :
    surface = surfArr[i]
    ax = fig.add_subplot( 231+i,projection='3d')
    ax.set(xlim=maxmin, ylim=maxmin, zlim=maxmin )
    ax.view_init(elev, azim)
    ax.set_title(str(i+1)+'. '+ titstrg[i], fontsize='x-small')
    ax.add_collection3d(surface)
    ax.set_axis_off()
fig.tight_layout()
plt.show()

