from matplotlib import pyplot as plt
import s3dlib.surface as s3d

# 1. Define functions to examine ....................................
# 2. Setup and map surfaces .........................................

earth = s3d.SphericalSurface(6)
earth.map_geom_from_image('data/elevation.png',0.1)
earth.shade()

# 3. Construct figure, add surfaces, plot ...........................

fig = plt.figure(figsize=plt.figaspect(1) )

ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1))
ax.view_init(elev=30, azim=120)

ax.add_collection3d(earth)

plt.show()