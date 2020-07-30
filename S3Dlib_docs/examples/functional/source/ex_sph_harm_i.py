import numpy as np
from scipy import special as sp
from matplotlib import pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Two methods of Representing f(θ, φ), imaginary

# 1. Define functions to examine ....................................

def sphHar(rtp) :
    r, theta, phi = rtp
    m, l = 2,3
    r = sp.sph_harm(m, l, theta, phi).imag
    return r, theta, phi

def sphHar_absR(rtp) :
    r, theta, phi = sphHar(rtp)
    return np.abs(r), theta, phi

# 2. Setup and map surfaces .........................................
rez = 5
binmap = cmu.binary_cmap('goldenrod','darkcyan',)

sph_23 = s3d.SphericalSurface(rez, basetype='octa', cmap='BrBG')
sph_23.map_cmap_from_op( lambda rtp : sphHar(rtp)[0] ).shade(0.8) 

sph_23_pos = s3d.SphericalSurface(rez, basetype='octa', cmap=binmap)
sph_23_pos.map_cmap_from_op( lambda rtp : sphHar(rtp)[0] )
sph_23_pos.map_geom_from_op(sphHar_absR).shade()

# 3. Construct figure, add surfaces, and plot .....................

fig = plt.figure(figsize=plt.figaspect(0.5/1.2))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')
ax1.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1) )
ax2.set(xlim=(-.3,.3), ylim=(-.3,.3), zlim=(-.3,.3) )
ax1.xaxis.set_major_locator(LinearLocator(5))
ax1.yaxis.set_major_locator(LinearLocator(5))
ax1.zaxis.set_major_locator(LinearLocator(5))
ax2.xaxis.set_major_locator(LinearLocator(5))
ax2.yaxis.set_major_locator(LinearLocator(5))
ax2.zaxis.set_major_locator(LinearLocator(5))
ax1.set_title('Color representation\n'+r'rgb=f($\theta$,$\varphi$), r=1')
ax2.set_title('Geometric representation\n'+r'R=|f($\theta$,$\varphi$)|')
plt.colorbar(sph_23.cBar_ScalarMappable, ax=ax1,  shrink=0.6 )
plt.colorbar(sph_23_pos.cBar_ScalarMappable, ax=ax2,  shrink=0.6 )

ax1.add_collection3d(sph_23)
ax2.add_collection3d(sph_23_pos)

fig.tight_layout()
plt.show()