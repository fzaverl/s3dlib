import numpy as np
import matplotlib.pyplot as plt
import s3dlib.surface as s3d

#.. Imaginary Earth

# 1. Define functions to examine ....................................

real = True
imaginary =  not real

def sqrt_Z(rtz, isReal) :
    r,t,z = rtz
    T=2*t
    if isReal :  Z = np.sqrt(r)*np.cos(T/2)
    else :       Z = np.sqrt(r)*np.sin(T/2)
    return r,T,Z

# 2. Setup and map surfaces .........................................

surface_1 = s3d.PolarSurface(6)
surface_1.map_color_from_image("data/earth.png")
surface_1.transform(s3d.eulerRot(115,0))
surface_1.map_geom_from_op( lambda rtz : sqrt_Z(rtz,imaginary) ).shade(.2,direction=[1,1,1])

# 3. Construct figure, add surfaces, and plot .....................

minmax = (-.8,.8) 
fig = plt.figure(figsize=plt.figaspect(1), facecolor='black')
ax1 = fig.add_subplot(111, projection='3d')
ax1.view_init(20, 205)
ax1.set(xlim=minmax, ylim=minmax, zlim=minmax)
ax1.set_facecolor('black')
ax1.set_title('Imaginary Earth',color='white')
ax1.set_axis_off()

ax1.add_collection3d(surface_1)

fig.tight_layout()
plt.show()