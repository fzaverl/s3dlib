from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Different Sub-surface Type

# 1 & 2. Define functions, setup and map surfaces ...................
rez=7

exterior = s3d.SphericalSurface(rez, basetype='octa')
exterior.map_color_from_image('data/blue_marble.png')
exterior.clip( lambda xyz : xyz[0]<0 , usexyz=True )

interior = s3d.PolarSurface(5, basetype='squ', cmap="hot")
interior.map_cmap_from_op( lambda rtz : 1-rtz[0]  )
interior.transform( [ [0,0,-1], [0,1,0], [1,0,0] ] )

surface = interior + exterior

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1), facecolor='black' )
desc = str(surface) + '\n' + str(exterior) + '\n' + str(interior)
fig.text(0.975,0.975, desc, ha='right', va='top', 
        fontsize='smaller', multialignment='right', color='white')
ax = fig.add_subplot(111, projection='3d')
minmax = (-0.75,0.75)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()
ax.set_facecolor('black')

ax.add_collection3d(surface)

fig.tight_layout()
plt.show()