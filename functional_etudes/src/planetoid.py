import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. #1 - influenced by M.C.Escher - Double Planetoid
#   https://mcescher.com/gallery/mathematical/

# 1. Define function to examine .....................................
grcolor = [0.890, 0.855, 0.824]
olcolor = [0.780, 0.753, 0.643]
bucolor = [0.239, 0.255, 0.298]
ttlcolor =  [0.502, 0.200, 0.278, 0.1] 
fgbgcolor = [0.886, 0.835, 0.773]
elev, azim = 58,25
light_direction = [1,-1,1 ]
illum = s3d.rtv(light_direction,elev,azim)

np.random.seed(2)

def randfunc(xyz) :
    r,t,z = s3d.SphericalSurface.coor_convert(xyz,False)
    sigma = 0.2*np.power(r,3)
    R = r + sigma*np.random.rand( len(r) )
    XYZ = s3d.SphericalSurface.coor_convert([R,t,z],True)
    return XYZ

def Tetrahedron(rez,color) :
    t = s3d.SphericalSurface(basetype='tetra')
    surface = s3d.Surface3DCollection(t.vertexCoor, t.fvIndices, facecolor=color)
    surface.triangulate(rez)
    surface.name = 'tetrahedron'
    surface.map_geom_from_op(randfunc)
    surface.set_facecolor(color)
    return surface

# 2. Setup and map surfaces .........................................
rez = 6

surface1 = Tetrahedron(rez,olcolor)
surface2 = Tetrahedron(rez,grcolor)
surface2.transform(rotate=s3d.eulerRot(180,180))
surface = surface1 + surface2
surface.shade(direction=illum)

disk = s3d.PolarSurface(4,color=bucolor)
disk.transform(rotate=s3d.eulerRot(-elev,-azim),scale=1.4)

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=plt.figaspect(1), facecolor=fgbgcolor)
text = fig.text(0.97, 0.87, 'S3Dlib.org', color=ttlcolor, ha='right',
    rotation=90, fontsize=55, fontweight='bold'  )
ax = plt.axes(projection='3d')
maxmin = ( -0.9,0.9 )
ax.set(xlim=maxmin, ylim=maxmin, zlim=maxmin )
ax.set_axis_off()
ax.view_init(elev, azim)
ax.set_facecolor(fgbgcolor)
ax.set_proj_type('ortho')

ax.add_collection3d(surface)
ax.add_collection3d(disk)

fig.tight_layout()
plt.show()