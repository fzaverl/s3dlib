from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Guides: Base Surfaces - Split Geometries

# 1. Define functions to examine ....................................

delta_z, delta_r, delta_t = -0.2, -0.2, -0.3

def planarSplit(xyz) :
    x,y,z = xyz
    d = 0.9
    v = x*y
    q = np.abs(v)
    s = np.sign(v)
    z = d*s*np.power(q,1/5)
    return x,y,z

def polarSplit(rtz) :
    r,t,z = rtz
    Z = (t/np.pi - 1)*delta_z/2
    return r,t,Z

def cylinSplit(rtz) :
    r,t,z = rtz
    T = t + (t/np.pi - 1)*delta_t/2
    R = r + (t/np.pi - 1)*delta_r/2
    return R,T,z

def sphereSplit(rtp) :
    r,t,p = rtp
    T = t + (t/np.pi - 1)*delta_t/2
    return r,T,p

basetype = OrderedDict()
basetype['polar'] =       [ s3d.PolarSurface,        polarSplit, 'squ_s', 'hex_s', 'squ_c', 'hex_c']
basetype['cylindrical'] = [ s3d.CylindricalSurface,  cylinSplit, 'tri_s', 'squ_s']
basetype['spherical'] =   [ s3d.SphericalSurface,    sphereSplit,'octa_s', 'cube_s','octa_c', 'cube_c']
# of subplots:       1        2           3           4           5           6
figlay =  [ [0,0], [1,1],   [1,2],      [1,3],      [2,2],      [2,3],      [2,3] ]
figsize = [ [0,0], [1,1], [4.8, 2.3], [7.2, 2.3], [4.8, 4.8], [7.2, 4.8], [7.2, 4.8] ]

rez=3
crdb = cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75])
minmax, ticks = (-1,1) , [-1,0,1]

for surface_type, surface_bases in basetype.items() :
    numbtypes = len(surface_bases) -2
    btypes = surface_bases[2:]
    i = 1
    rows, colm = figlay[numbtypes][0], figlay[numbtypes][1]
    width, height = figsize[numbtypes][0], figsize[numbtypes][1]
    fig = plt.figure(figsize=(width,height))
    for typeString in btypes :

        # Setup and map surfaces .........................................
        surObj = surface_bases[0](rez,basetype=typeString, linewidth=0.2)
        surObj.map_geom_from_op(surface_bases[1])
        surObj.map_cmap_from_normals(cmap=crdb, direction=[1,1,1]).shade(0.5)
        surObj.set_edgecolor( [0.25,0.15,0] )

        # Construct figures, add surfaces, plot ......................
        ax = fig.add_subplot(rows, colm,i, projection='3d')
        ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
        ax.set_proj_type('ortho')
        ax.set_xlabel('X axis', labelpad = 0)
        ax.set_ylabel('Y axis', labelpad = 0)
        ax.set_zlabel('Z axis', labelpad = 0)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_zticks(ticks)
        i += 1
        ax.text(0,1,1,'  '+typeString, color='tab:red', horizontalalignment='left', verticalalignment='bottom')
        ax.view_init(elev=45, azim=-25)

        ax.add_collection3d(surObj)

plt.show()

