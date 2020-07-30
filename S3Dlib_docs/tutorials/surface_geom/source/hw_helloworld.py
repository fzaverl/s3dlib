import matplotlib.pyplot as plt
import s3dlib.surface as s3d

'''
HELLO WORLD
'''

# Setup surface ................................................

surface = s3d.SphericalSurface()
surface.shade()

# Construct figure, add surface, plot ..........................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
ax.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1))

ax.add_collection3d(surface)

plt.show()