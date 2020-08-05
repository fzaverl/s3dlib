import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

#.. Guides: Vector Fields -
#   surface with vector plots

# setup functions ===================================================

delta_1, delta_2 = 0,0

def setdelta(rez, size_1, size_2) :
    global delta_1
    global delta_2
    sigma = 0.4
    n_verts = s3d.SphericalSurface.fev(rez,basetype='octa')[2]
    delta_1 = np.full( n_verts, size_1) + sigma*np.random.rand( n_verts)
    delta_2 = np.full( n_verts, size_2) + sigma*np.random.rand( n_verts)
    return

def original(rtp) :
    r,t,p = rtp
    r = delta_1
    return r,t,p

def shifted(rtp) :
    r,t,p = rtp
    r = delta_2
    return r,t,p

def vect_field(rtp) :
    r,t,p = rtp
    val = 0.2
    u = np.full(len(r),val)
    v = np.full(len(r),val)
    w = np.full(len(r),val)
    return u,v,w

# setup surface =====================================================

def set_normals(ax) :
    rez = 1
    setdelta(rez, 0.8, 1)
    surface = s3d.SphericalSurface(rez,basetype='octa', color=[1,0,0,.2])
    surface.map_geom_from_op(original).shade(direction=[-1,1,1])
    vf = surface.facenormals( scale=0.4 )

    ax.add_collection3d(surface)
    ax.add_collection3d(vf)
    return

def set_vectorfield(ax) :
    rez = 2
    setdelta(rez, 0.7, 1)
    surface = s3d.SphericalSurface(rez,basetype='octa', color=[1,0,0,.2])
    surface.map_geom_from_op(original)
    surface.shade(direction=[-1,1,1])
    vf = surface.vectorfield_from_op(vect_field, scale=1)

    ax.add_collection3d(surface)
    ax.add_collection3d(vf)
    return

def set_dispfield(ax) :
    rez = 1
    setdelta(rez, 0.7, 1)
    surface = s3d.SphericalSurface(rez,basetype='octa', color=[1,0,0,.2])
    surface.map_geom_from_op(original)
    surface.shade(direction=[-1,1,1])
    vf = surface.dispfield_from_op(shifted,returnxyz=False)
    surface2 = s3d.SphericalSurface(rez,basetype='octa', edgecolor=[0,0,1,.2])
    surface2.map_geom_from_op( shifted )
    edges = surface2.edges
    edges.set_color([0,0,1,.2])

    ax.add_collection3d(surface)
    ax.add_collection3d(vf)
    ax.add_collection3d(edges)
    return

def set_surface(ax) :
    rez = 1
    setdelta(rez, .6, .4)
    surface = s3d.SphericalSurface(rez,basetype='octa', color=[1,0,0,.2])
    surface.map_geom_from_op(original).shade(direction=[-1,1,1])
    surface.transform(translate=(0.5,-.5,0.8))

    surface2 = s3d.SphericalSurface(rez,basetype='octa', color=[0,0,1,.2])
    surface2.map_geom_from_op( shifted ).shade(direction=[-1,1,1])
    surface2.transform(translate=(.5,1.2,-.3))

    vf = surface.vectorfield_to_surface(surface2,scale=0.3,alr=.2)
    vf2 = surface.vectorfield_to_surface(surface2,color=[0,0,1,.2],alr=0)

    ax.add_collection3d(surface)
    ax.add_collection3d(surface2)
    ax.add_collection3d(vf)
    ax.add_collection3d(vf2)
    return

# construct fig, plot ===============================================
funcArr = [set_normals, set_vectorfield, set_dispfield, set_surface]
minmax = (-1,1)

for i in range(4) :
    fig = plt.figure(figsize=plt.figaspect(1))
    ax = plt.axes(projection='3d')
    ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1))
    ax.set_axis_off()
    s3d.setupAxis(ax, length=1.75,negaxis=False)
    ax.view_init(30, 40)
    funcArr[i](ax)
    fig.tight_layout()

plt.show()