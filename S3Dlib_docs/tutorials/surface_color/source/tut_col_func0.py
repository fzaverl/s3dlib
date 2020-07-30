import matplotlib.pyplot as plt
import numpy as np

# 1. Define function to examine ....................................

def f_func(x,y) :
    r = np.sqrt( x**2 + y**2)
    Z = np.sin( 6.0*r )/2
    return Z
    
# 2. Setup and map plane ...........................................

x = np.linspace(-1,1.0,64)
y = np.linspace(-1,1.0,64)
xx, yy = np.meshgrid(x, y)
z = f_func(xx,yy)

# 3. Construct figure, add surface, and plot ......................

fig = plt.figure(figsize=(2.65,2.65))
ax = plt.axes()

ax.imshow(z, cmap='RdYlGn')

plt.show()
