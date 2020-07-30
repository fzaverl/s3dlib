from matplotlib import pyplot as plt
import numpy as np
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Texture Surface with Geometric Mapping

# 1. Define functions to examine ....................................

def randfunc(rtz) :
    r,t,z = rtz
    sigma = 0.01
    R = r + sigma*np.random.rand( len(r) )
    return R,t,z

def ava_func(rtp) :
    r,t,p = rtp
    delta, tmax = 0.1, 0.7
    P = 2*p/tmax
    Rt = (1-delta) + delta*np.cos(P)
    R = r*np.where(P>2*np.pi, 1.0 , Rt)
    return R,t,p

def ava_edge(rtz) :
    r,t,z = rtz
    tb = 2*np.pi - t
    T = np.where(t<np.pi,t,tb)
    # use the ava_func only to calc R
    rr = np.full(len(r),0.98)
    R = r*ava_func([rr, z, T])[0]
    return R,t,z

# 2. Setup and map surfaces .........................................
rez=6
cm_skin = cmu.rgb_cmap_gradient('darkolivegreen', 'olive' )
cm_seed = cmu.rgb_cmap_gradient('saddlebrown', 'sienna' )
cm_intr = cmu.rgb_cmap_gradient('saddlebrown', 'darkkhaki' )

skin = s3d.SphericalSurface(rez,basetype='octa')
skin.map_geom_from_op(randfunc)
skin.map_cmap_from_op( lambda rtz : rtz[0], cm_skin )
skin.map_geom_from_op( lambda rtp : ava_func(rtp) )
skin.shade().hilite(.3,direction=[0,1,1])
skin.clip( lambda xyz : xyz[0]>0 , usexyz=True )

interior = s3d.PolarSurface(rez, basetype='squ')
interior.map_cmap_from_op( lambda rtz : rtz[0] , cm_intr  )
interior.map_geom_from_op(ava_edge)
interior.transform( [ [0,0,1], [0,1,0], [1,0,0] ] )

seed =  s3d.SphericalSurface(rez,basetype='octa')
seed.map_cmap_from_op( lambda rtp : randfunc(rtp)[0], cm_seed )
seed.transform(scale=.4,translate=[0,0,-.1]).shade()

avacado = interior + skin + seed

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1))
info = str(avacado) 
fig.text(0.975,0.975,info, ha='right', va='top', fontsize='smaller', multialignment='right')
ax = fig.add_subplot(111, projection='3d')
maxmin = (-.8,.8)
ax.set(xlim=maxmin, ylim=maxmin, zlim=maxmin )
ax.set_axis_off()
ax.view_init(elev=10, azim=125)

ax.add_collection3d(avacado)

fig.tight_layout()
plt.show()