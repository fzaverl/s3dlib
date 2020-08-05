import matplotlib.pyplot as plt
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Guides: Orientation - Object Rotations
#   euler_rotations

rez = 5
theta, phi = 60, 45
initColor, testColor = [0.4,0.4,0.8, 0.2] , [0.4,.7,0.4, 0.2]

Xcon = True
comboSurface, title, lineColl = [], [], []
for i in range(2) :
    r = s3d.PolarSurface(rez, color=initColor)
    r.transform(rotate=s3d.eulerRot(theta,0,useXconv=Xcon))
    s = s3d.PolarSurface(rez, color=testColor)
    s.transform(rotate=s3d.eulerRot(theta,phi,useXconv=Xcon))
    comboSurface.append( r + s )
    lineColl.append(s.get_transformAxis(1.25))
    t = r'$\theta$ = {}, $\phi$ = {} '.format(theta,phi)
    t2 = '\nX-convention (zxz)'
    if not Xcon : t2 = '\nY-convention (zyz)'
    title.append(t+t2)
    Xcon = not Xcon

# 3. Construct figure, add dasurfaces, and plot ....................
fig = plt.figure(figsize=(7,3.5))
for i in range(2) :
    ax = fig.add_subplot(1,2,i+1, projection='3d')
    ax.set_title(title[i])
    ax.set_axis_off()
    ax.add_collection3d(comboSurface[i])
    ax.add_collection3d(lineColl[i])
    s3d.standardAxis( ax )
fig.tight_layout()
plt.show()
