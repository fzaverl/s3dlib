import numpy as np
from scipy import special as sp
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

# 1. Define functions to examine ....................................
def sphHar(rtp, m, n) :
    r, theta, phi = rtp
    R = sp.sph_harm(m, n, theta, phi).real
    return R, theta, phi

def sphHar_absR(rtp, m, n) :
    r, theta, phi = sphHar(rtp, m, n)
    return np.abs(r), theta, phi

# 2 & 3. Construct figure & axes, add surfaces, show ................

fig = plt.figure(figsize=plt.figaspect(1))
for m in range(3) :
    for n in range(1,4) :    
        i = 3*m + n
        if i==7 : continue
        ax = fig.add_subplot(3,3,i, projection='3d')
        ax.set_title('('+str(m)+','+str(n)+')' , fontsize='large')
        
        surface = s3d.SphericalSurface(5, cmap="RdYlGn")
        surface.map_cmap_from_op( lambda rtp : sphHar(rtp,m,n)[0] )
        surface.map_geom_from_op( lambda rtp : sphHar_absR(rtp,m,n) )
        s3d.auto_scale(ax,surface,rscale=0.6).set_axis_off()
        
        ax.add_collection3d(surface.shade())

fig.tight_layout()
plt.show()