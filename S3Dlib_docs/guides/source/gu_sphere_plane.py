import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. sphere vs planar, polar, cylindrical-sphere vertices

# 1. Define function to examine .....................................

def plane2sphere(xyz) :
    x,y,z = xyz
    T = -np.pi*(x+1)
    P = np.pi*(y+1)/2
    R = np.ones(len(z)).astype(float)
    XYZ = s3d.SphericalSurface.coor_convert([R,T,P],True)
    return XYZ

def disk2sphere(rtz) :
    r,t,z = rtz
    P = np.pi*r
    R = np.ones(len(z)).astype(float)
    XYZ = s3d.SphericalSurface.coor_convert([R,t,P],True)
    return XYZ

def cylinder2sphere(rtz) :
    r,t,z = rtz
    P = np.pi*(z+1)/2
    t = -t
    XYZ = s3d.SphericalSurface.coor_convert([r,t,P],True)
    return XYZ

# 2. Setup and map surfaces .........................................
rez = 3
c = 'khaki'

plane = s3d.PlanarSurface(rez+1,basetype='oct1',color=c)
x,y,z = plane.vertices
plane.map_geom_from_op( plane2sphere ).shade().hilite(.5)

disk = s3d.PolarSurface(rez+1,basetype='hex',color=c)
rp,tp,zp = s3d.PolarSurface.coor_convert(disk.vertices)
disk.map_geom_from_op( disk2sphere,True ).shade().hilite(.5)

cylinder = s3d.CylindricalSurface(rez,basetype='tri',color=c)
rc,tc,zc = s3d.CylindricalSurface.coor_convert(cylinder.vertices)
cylinder.map_geom_from_op( cylinder2sphere,True ).shade().hilite(.5)

sphere = s3d.SphericalSurface(rez,basetype='icosa',color=c).shade().hilite(.5)
rs,ts,ps = s3d.SphericalSurface.coor_convert(sphere.vertices)

surfaces = [
    [ plane,    [x,y],    '(X,Y)',                'Planar' ],
    [ disk,     [tp,rp], r'($\theta$,R)',         'Polar' ],
    [ cylinder, [tc,zc], r'($\theta$,Z)',         'Cylindrical' ],
    [ sphere,   [ts,ps], r'($\theta$,$\varphi$)', 'Spherical' ] ]

# 3. Construct figure, add surface, plot ............................

# Figure 1 Spherical surface views .........................
fig = plt.figure(figsize=(6,6))
maxmin = ( -0.8,0.8 )
for i in range(len(surfaces)) :
    surface = surfaces[i]
    ax = fig.add_subplot(221+i, projection='3d')
    ax.set(xlim=maxmin, ylim=maxmin, zlim=maxmin )
    title = surface[3] + 'Surface, ' + str(len(surface[0].vertices[0]))
    ax.set_title( title )
    ax.set_axis_off()
    ax.set_proj_type('ortho')
    ax.view_init(0)
    ax.add_collection3d(surface[0])

fig.tight_layout()

# Figure 2 - Vertex distribution plots .....................
fig = plt.figure(figsize=(6,6))
for i in range(len(surfaces)) :
    surface = surfaces[i]
    ax = fig.add_subplot(221+i)
    title = surface[2] +' '+ surface[3] +' Vertices'
    ax.set_title( title , fontsize='medium' )
    x,y = surface[1]
    ax.scatter(x,y,s=1)

# ..........................................................
plt.show()