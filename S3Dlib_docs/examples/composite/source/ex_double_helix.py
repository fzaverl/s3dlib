import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import LinearLocator
from matplotlib.colors import ListedColormap 
import s3dlib.surface as s3d
import s3dlib.cmap_utilities as cmu

#.. Angular 4-Color Color Map

# 1. Define functions to examine ....................................

code = "AGATTCAGTACTAGAGTCCGAGATAGGAT"
k = 3                 # number of twists
p1, p2 = 13.0, 21.0   # double helix spacing
delta = np.pi*p1/(p1+p2)

def basePairs( codeString ) :
    innerMap = { 'A':0, 'G':1, 'T':2, 'C':3 }
    outerMap = { 'T':0, 'C':1, 'A':2, 'G':3 }
    posMap, negMap = [], []
    for code in codeString : 
        posMap.append(outerMap[code])
        negMap.append(innerMap[code])
    return np.array( [posMap,negMap] )

indexMap = basePairs(code)

def getColorIndex(rtz) :
    r,t,z = rtz
    codeLen = len(indexMap[1])
    posNeg = np.where(z>0,np.zeros(len(z)),np.ones(len(z)))
    codeMap = indexMap[posNeg.astype(int)]
    tNorm = (codeLen*(t/(2*np.pi))).astype(int)
    tColor = np.choose( tNorm, np.transpose(codeMap) )  
    return tColor

def nacidShape(rtz) :
    r,t,z = rtz
    Rmin = np.cos(delta)
    Z = k*t - k*np.pi
    T = k*t + z*delta   
    R = Rmin/np.cos(z*delta)
    return R,T,Z

def riboShape(rtz, isPos=True) :
    r,t,z = rtz
    ratio = .1
    Z = ratio*np.sin(z*np.pi)
    R = r + ratio*np.cos(z*np.pi)
    T = k*t
    Z = k*Z + T - k*np.pi
    ofset = -delta
    if isPos : ofset = delta
    return R,T+ofset,Z
    
# 2. Setup and map surfaces .........................................
rez= 6
cmap = ListedColormap(['red','yellow','green','blue'])

nacid = s3d.CylindricalSurface(rez, basetype='squ_s')
nacid.map_cmap_from_op(getColorIndex,cmap)
nacid.map_geom_from_op( nacidShape ).shade(0.3)

riboA = s3d.CylindricalSurface(rez, basetype='squ_s', color='darkgoldenrod')
riboA.map_geom_from_op( lambda rtz : riboShape(rtz) ).shade().hilite(0.9)
riboB = s3d.CylindricalSurface(rez, basetype='squ_s', color='darkgoldenrod')
riboB.map_geom_from_op( lambda rtz : riboShape(rtz,False) ).shade().hilite(0.9)

dna = nacid + riboA + riboB

# 3. Construct figure, add surfaces, and plot ......................

fig = plt.figure(figsize=plt.figaspect(1), facecolor='black')
fig.text(0.975,0.975,str(dna), ha='right', va='top', fontsize='smaller', multialignment='right', color='white')
tL, tB, tD = 0.80, 0.28, 0.137
fig.text(tL,tB,     'Adenine',  color='white', ha='right')
fig.text(tL,tB+1*tD,'Guanine',  color='white', ha='right')
fig.text(tL,tB+2*tD,'Thymine',  color='white', ha='right')
fig.text(tL,tB+3*tD,'Cytosine', color='white', ha='right')
ax = plt.axes(projection='3d', facecolor='black')
ax.set(xlim=(-1.5,1.5), ylim=(-1.5,1.5), zlim=(-7,7) )
plt.colorbar(nacid.cBar_ScalarMappable, ax=ax, shrink=0.6, pad=0.1 )
ax.set_axis_off()
ax.view_init(elev=10)

ax.add_collection3d(dna)

fig.tight_layout()
plt.show()