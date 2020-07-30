import matplotlib.pyplot as plt
from collections import OrderedDict
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

basetype = OrderedDict()

basetype['planar'] =      [ s3d.PlanarSurface,      'quad','oct1','oct2']
basetype['polar'] =       [ s3d.PolarSurface,       'squ','hex']
basetype['cylindrical'] = [ s3d.CylindricalSurface, 'tri', 'tri2','squ','squ2']
basetype['spherical'] =   [ s3d.SphericalSurface,   'tetra','octa','icosa','cube','dodeca']

minmax = (-1,1)
cmu.rgb_cmap_gradient([0.25,0.15,0],[1,.9,.75],'cardboard')

for surface_type, surface_bases in basetype.items() :
    numbtypes = len(surface_bases) -1
    fig = plt.figure(figsize=plt.figaspect(1/numbtypes))
    btypes = surface_bases[1:]
    i = 1
    for typeString in btypes :
        # Setup and map surfaces .........................................
        
        surObj = surface_bases[0](basetype=typeString)       
        surObj.map_cmap_from_normals(cmap='cardboard')
        surObj.set_edgecolor( [0.25,0.15,0] )

        # Construct figures, add surfaces, plot ......................
        ax = fig.add_subplot(1,numbtypes,i, projection='3d')
        ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
        ax.set_title(typeString)
        ax.set_proj_type('ortho')
        i += 1
        
        ax.add_collection3d(surObj)

plt.show()

