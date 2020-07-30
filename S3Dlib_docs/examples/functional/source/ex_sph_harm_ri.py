import numpy as np
from scipy import special as sp
from matplotlib import pyplot as plt
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Two methods of Representing f(θ, φ), real & imaginary (abs values)

# 1. Define functions to examine ....................................
real = True
imaginary =  not real

def sphHar(rtp, isReal) :
    r, theta, phi = rtp
    m, l = 2,3
    sph = sp.sph_harm(m, l, theta, phi)
    if isReal : r = sph.real
    else :      r = sph.imag
    return r, theta, phi

def sphHar_absR(rtp, isReal) :
    r, theta, phi = sphHar(rtp, isReal)
    return np.abs(r), theta, phi

# 2. Setup and map surfaces .........................................
rez = 6
binmap = cmu.binary_cmap('goldenrod','darkcyan',)

sph_23_imag = s3d.SphericalSurface(rez, basetype='octa', cmap=binmap)
sph_23_imag.map_cmap_from_op( lambda rtp : sphHar(rtp,imaginary)[0] )
sph_23_imag.map_geom_from_op( lambda rtp : sphHar_absR(rtp,imaginary) ).shade()

sph_23_real = s3d.SphericalSurface(rez, basetype='octa', cmap=cmu.binary_cmap())
sph_23_real.map_cmap_from_op( lambda rtp : sphHar(rtp,real)[0] )
sph_23_real.map_geom_from_op( lambda rtp : sphHar_absR(rtp,real) ).shade()

sph_23_pos = sph_23_real + sph_23_imag

# 3. Construct figure, add surfaces, and plot .....................

fig = plt.figure(figsize=plt.figaspect(1))
ax = fig.add_subplot(111, projection='3d')
ax.set(xlim=(-.3,.3), ylim=(-.3,.3), zlim=(-.3,.3) )
ax.set_title('Real and Imaginary\n'+r'R=|f($\theta$,$\varphi$)|')
ax.set_axis_off()

ax.add_collection3d(sph_23_pos)

fig.tight_layout()
plt.show()