from collections import OrderedDict
import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Guides: Base Surfaces - General Geometries

basetype = OrderedDict()
basetype['planar'] =      [ s3d.PlanarSurface,      'quad','oct1','oct2']
basetype['polar'] =       [ s3d.PolarSurface,       'squ','hex']
basetype['cylindrical'] = [ s3d.CylindricalSurface, 'tri', 'tri2','squ','squ2']
basetype['spherical'] =   [ s3d.SphericalSurface,   'tetra','octa','icosa','cube','dodeca']
figlay = [ [0,0], [1,1], [1,2], [1,3], [2,2], [2,3], [2,3] ]
#                    1         2           3            4           5           6
figsize = [ [0,0], [1,1], [4.8, 2.3], [7.2, 2.3], [4.8, 4.8], [7.2, 4.8], [7.2, 4.8] ]

cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboard')
minmax, ticks = (-1,1) , [-1,0,1]

for surface_type, surface_bases in basetype.items() :
    numbtypes = len(surface_bases) -1
    btypes = surface_bases[1:]
    i = 1
    rows, colm = figlay[numbtypes][0], figlay[numbtypes][1]
    width, height = figsize[numbtypes][0], figsize[numbtypes][1]
    fig = plt.figure(figsize=(width,height))
    for typeString in btypes :
        
        # Setup and map surfaces .........................................
        surObj = surface_bases[0](basetype=typeString)       
        surObj.map_cmap_from_normals(cmap='cardboard', direction=[1,1,1])
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

        ax.add_collection3d(surObj)

plt.show()

