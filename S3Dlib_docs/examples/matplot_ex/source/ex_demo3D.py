import matplotlib.pyplot as plt
import s3dlib.surface as s3d

#.. Matplotlib Examples: 3D surface (solid color)

# Setup surface ................................................

surface = s3d.SphericalSurface(4).shade()
surface.transform(scale=10)
#surface.map_cmap_from_normals()

# Construct figure, add surface, plot ..........................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
ax.set(xlim=(-10,10), ylim=(-10,10), zlim=(-10,10))

ax.add_collection3d(surface)

plt.show()