import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. #1 - influenced by M.C.Escher - Tetrahedral Planetoid
#   https://mcescher.com/gallery/mathematical/

# 1. Define function to examine .....................................
ttlcolor =  [0.502, 0.200, 0.278, 0.1] 
mncolor =   [ 0.867, 0.816, 0.780 ]
gbcolor =   [ 0.224, 0.208, 0.204 ]
fgbgcolor = [ 0.863, 0.800, 0.737 ]
elev, azim = 35, 0
light_direction = [1,.3,1 ]
illum = s3d.rtv(light_direction,elev,azim)


def Tetrahedron(rez,color) :
    t = s3d.SphericalSurface(basetype='tetra')
    surface = s3d.Surface3DCollection(t.vertexCoor, t.fvIndices, color=color)
    surface.triangulate(rez)
    surface.name = 'tetrahedron'
    surface.set_color(color)
    return surface

outerdpth, innerdpth = 0.48,0.4

def pushIn(xyz) :
    r,t,p = s3d.SphericalSurface.coor_convert(xyz)
    minR = np.full(len(r),outerdpth)
    delta = np.power(  np.abs(minR-r) , 0.65 )
    R = np.where(r<minR,r - delta, r)
    XYZ = s3d.SphericalSurface.coor_convert([R,t,p],True)
    return XYZ

def pushOut(xyz) :
    r,t,p = s3d.SphericalSurface.coor_convert(xyz)
    minR = np.full(len(r),innerdpth)
    R = np.where(r<minR,minR,r)
    XYZ = s3d.SphericalSurface.coor_convert([R,t,p],True)
    return XYZ

def randfunc(xyz) :
    r,t,p = s3d.SphericalSurface.coor_convert(xyz)
    sigMx = 0.3
    a = sigMx/(1.0-outerdpth)
    sigma = a*( r - np.full(len(r),outerdpth) )
    R = np.where(r>outerdpth, r + sigma*np.random.rand( len(r) ), r)
    XYZ = s3d.SphericalSurface.coor_convert([R,t,p],True)
    return XYZ

# 2. Setup and map surfaces .........................................
rez, srez = 6, 4
cmap = cmu.hsv_cmap_gradient([0.16,.3,.3,1],[0.16,.3,.3,.0])

surface = Tetrahedron(rez,mncolor)
surface.map_geom_from_op(pushIn)
surface.map_geom_from_op(randfunc)
surface.map_geom_from_op(pushOut)
surface.shade(direction=illum)

sphere = s3d.SphericalSurface(srez).transform(scale=.7)
normillum = s3d.elev_azim_2vector(elev,azim)
sphere.map_cmap_from_normals(cmap,direction=normillum)
surface = surface + sphere

# 3. Construct figure, add surface, plot ............................

fig = plt.figure(figsize=plt.figaspect(1), facecolor=fgbgcolor)
text = fig.text(0.97, 0.87, 'S3Dlib.org', color=ttlcolor, ha='right',
    rotation=90, fontsize=55, fontweight='bold'  )
ax = plt.axes(projection='3d')
ax.set_proj_type('ortho')
minmax = (-.6,.6)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_facecolor(gbcolor)
ax.set_axis_off()
ax.view_init(elev,azim)

ax.add_collection3d(surface)

fig.tight_layout(rect=(.05,.05,.95,.95))
plt.show()