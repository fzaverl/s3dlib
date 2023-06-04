from matplotlib import pyplot as plt
import s3dlib.surface as s3d

# 2. Setup and map surfaces .........................................
earth = s3d.SphericalSurface(6)
earth.map_color_from_image('data/earth.png')
earth.map_geom_from_image('data/elevation.png',0.06)
earth.transform(rotate=s3d.eulerRot(175,0) )

# 3. Construct figure, add surfaces, and plot ......................
fig = plt.figure(figsize=plt.figaspect(1), facecolor='k' )
ax = plt.axes(projection='3d', facecolor='k')
minmax = (-0.65,0.65)
ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax.set_axis_off()

ax.add_collection3d(earth.shade( contrast=1.7, direction=[1,0.8,1] ))

plt.show()