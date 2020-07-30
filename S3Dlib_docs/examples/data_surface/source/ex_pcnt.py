import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import s3dlib.surface as s3d

#.. Percentile Visualization

# 1. Define data to examine .........................................
np.random.seed(1)

def get_correlated_dataset(n, dependency, mu, scale):
    latent = np.random.randn(n, 3)
    dependent = latent.dot(dependency)
    scaled = dependent * scale       
    scaled_with_offset = scaled + mu
    # return x y z of the new, correlated dataset
    return scaled_with_offset[:, 0], scaled_with_offset[:, 1], scaled_with_offset[:, 2]

corr = np.array([ [0.85, -.15, 0.4], [-0.35, -0.65, 0.7], [-.4, 0.6, 1.0] ])
mu = 0,0,0
sigma = 1.35, 0.56 , 0.68
N = 400
x,y,z = get_correlated_dataset(N, corr, mu, sigma)
data = np.transpose([ x,y,z ])

# 2. Setup and map surfaces .........................................
prct = 0.95

surface = s3d.SphericalSurface(3, color=[0,0,0,0.05], linewidth=.5  )
disArr_a,t = surface.svd(data, prct)
surface.shade()
info = str(N) +', '+"{:.0%}".format(prct) + ', ' + '{:04.2f}'.format( t[0] )

colors_a = []
colorMap = cm.get_cmap('viridis')
for val in disArr_a :
    if val > 1 : colors_a.append([0.8,0,0])
    else : colors_a.append(colorMap(val))

# 3. Construct figures, add surfaces, and plot .......................

fig = plt.figure(figsize=plt.figaspect(1))
ax = plt.axes(projection='3d')
ax.set(xlim=[-3,3], ylim=[-3,3], zlim=[-3,3] )
ax.set_title(info, horizontalalignment='left')
ax.scatter(x,y,z, c=colors_a, edgecolor='k')

ax.add_collection3d(surface.get_transformAxis())
ax.add_collection(surface)

fig.tight_layout()
plt.show()