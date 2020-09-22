import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. #1 - influenced by M.C.Escher - Spirals
#   https://mcescher.com/gallery/mathematical/

# 1. Define function to examine .....................................
scolor =    [0.859, 0.788, 0.729]
gbcolor =   [0.506, 0.482, 0.435]
ttlcolor =  [0.502, 0.200, 0.278, 0.1] 
fgbgcolor = [0.867, 0.800, 0.729]
elev, azim = 45,30
light_direction = [0,-1,1 ]
illum = s3d.rtv(light_direction,elev,azim)
negillum = -1.0*np.array(illum)

def twisted_torus(rtz,rotate,isInner) :
    # modified function from documentation example 
    twists, width, rMax, rdMin, stretch = 5,.125,.7,.05, 1.4
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
    t1.map_geom_from_op( lambda rtz : twisted_torus(rtz,0.00,inner) )
    t2.map_geom_from_op( lambda rtz : twisted_torus(rtz,0.25,inner) )
    t3.map_geom_from_op( lambda rtz : twisted_torus(rtz,0.50,inner) )
    t4.map_geom_from_op( lambda rtz : twisted_torus(rtz,0.75,inner) )
    surface = t1 + t2 + t3 + t4
    return surface
# 2. Setup and map surfaces .........................................
rez = 8

surface1 = Spiralset(rez)
surface1.shade(direction=illum).hilite(.5,direction=illum)
surface2 = Spiralset(rez,True)
surface2.shade(direction=negillum).hilite(.5,direction=negillum)

surface = surface1 + surface2

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=(5,4) , facecolor=fgbgcolor)
text = fig.text(0.97, 0.87, 'S3Dlib.org', color=ttlcolor, ha='right',
    rotation=90, fontsize=45, fontweight='bold'  )
fig.text(0.11,0.12,str(surface),color=scolor ,fontsize='x-small')
ax = plt.axes(projection='3d')
maxmin = ( -0.9,0.9 )
ax.set(xlim=maxmin, ylim=(-.63,1.17), zlim=maxmin )
ax.set_axis_off()
ax.view_init(elev, azim)
ax.set_facecolor(gbcolor)

ax.add_collection3d(surface)

fig.tight_layout(rect=(.07,.07,.93,.93))
plt.show()