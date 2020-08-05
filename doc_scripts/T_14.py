import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Tutorial: Coordinates and Functions - Function Coordinates

arc_seg = 25
theta = np.pi/3
phi = np.pi/3.5

thetaArc = np.linspace(0.0,theta, arc_seg)
arcz = np.linspace(0, 0, arc_seg)
r = arcz + 1.0
arcx = r * np.cos(thetaArc)
arcy = r * np.sin(thetaArc)

axisLength = 1.5
fcolor0 = [0.7,0.7,0.6,0.4]
fcolor =  [0.7,0.7,0.6,0.2]
shade = .3
zoffset = 0.7
auxColor = 'silver'

def sphere(ax) :
    surface = s3d.SphericalSurface(4,facecolor=fcolor)
    surface.shade(depth=shade, direction=[-1,1,1])
    surface.set_edgecolor([0,0,0,0])

    phi_arc = np.linspace(0.0, phi, arc_seg)

    phi_z = np.cos(phi_arc)
    phi_x = np.sin(phi_arc)*np.cos(theta)
    phi_y = np.sin(phi_arc)*np.sin(theta)

    z = np.cos(phi)
    x = np.sin(phi)*np.cos(theta)
    y = np.sin(phi)*np.sin(theta)

    rfct = 0.4
    r_x = rfct*x
    r_y = rfct*y
    r_z = 0.65*z

    tfct = 1.3
    t_x = tfct*np.cos(theta/2)
    t_y = tfct*np.sin(theta/2)

    pfct = 0.3
    p_x = pfct*x
    p_y = pfct*y
    p_z = .85

    pfct = 0.85
    p_z = pfct*np.cos(phi/4)
    p_x = pfct*np.sin(phi/4)*np.cos(theta)
    p_y = pfct*np.sin(phi/4)*np.sin(theta)

    xcon,ycon = axisLength*np.cos(theta) , axisLength*np.sin(theta)

    arc_line = np.transpose( [ [0,0,0], [x,y,0] ] )
    ver_line = np.transpose( [ [x,y,0], [x,y,z] ] )
    rad_line = np.transpose( [ [0,0,0], [x,y,z] ] )
    arc_line_co_0 = np.transpose( [[x,y,0] , [np.cos(theta),np.sin(theta),0]] )
    arc_line_cont = np.transpose( [[np.cos(theta),np.sin(theta),0] , [xcon,ycon,0]] )
    ax.plot(arc_line[0], arc_line[1], arc_line[2], c='grey', linestyle=':')
    ax.plot(arc_line_co_0[0], arc_line_co_0[1], arc_line_co_0[2], c='grey', linestyle='-')
    ax.plot(arc_line_cont[0], arc_line_cont[1], arc_line_cont[2], c='black', linestyle='-')
    ax.plot(ver_line[0], ver_line[1], ver_line[2], c='grey', linestyle=':')
    ax.plot(rad_line[0], rad_line[1], rad_line[2], c='red', linestyle='-')
    ax.plot(arcx, arcy, arcz, c='blue', linewidth=1, linestyle='-')
    ax.plot(phi_x, phi_y, phi_z, c='green', linewidth=1, linestyle='-')

    ax.scatter(x,y,z,c='tab:red',s=40)
    ax.add_collection3d(surface)

    ax.text(t_x,t_y,0,r'$\theta$', fontweight='bold', fontsize='large', color = 'blue', horizontalalignment='center', verticalalignment='center')
    ax.text(p_x,p_y,p_z,r'$\phi$', fontsize='large', color = 'green', horizontalalignment='center', verticalalignment='center')
    ax.text(r_x,r_y,r_z,'r', fontsize='large', color = 'red', horizontalalignment='center', verticalalignment='center')

    return

def cylinder(ax) :
    surface = s3d.CylindricalSurface(4,facecolor=fcolor)
    surface.shade(depth=shade, direction=[-1,1,1])
    surface.set_edgecolor([0,0,0,0])
    
    x = np.cos(theta)
    y = np.sin(theta)
    z = zoffset

    xcon,ycon = axisLength*x , axisLength*y
    
    arc_line = np.transpose( [ [0,0,0], [x,y,0] ] )
    ver_line = np.transpose( [ [x,y,0], [x,y,z] ] )
    arc_line_cont = np.transpose( [[x,y,0] , [xcon,ycon,0]] )
    ax.plot(arc_line[0], arc_line[1], arc_line[2], c='red', linestyle='-')
    ax.plot(arc_line_cont[0], arc_line_cont[1], arc_line_cont[2], c='black', linestyle='-')
    ax.plot(ver_line[0], ver_line[1], ver_line[2], c='green', linestyle='-')
    ax.plot(arcx, arcy, arcz, c='blue', linewidth=1, linestyle='-')

    ax.scatter(x,y,z,c='tab:red',s=40)
    ax.add_collection3d(surface)

    ax.text(1.2,0.6,0,r'$\theta$', fontweight='bold', fontsize='large', color = 'blue', horizontalalignment='center', verticalalignment='center')
    ax.text(.7,.83,.55,'z', fontsize='large', color = 'green', horizontalalignment='center', verticalalignment='center')
    ax.text(.45,.2,0,'r', fontsize='large', color = 'red', horizontalalignment='center', verticalalignment='center')
    
    return

def disk(ax) :
    surface = s3d.PolarSurface(4,facecolor=fcolor0)
    surface.set_edgecolor([0,0,0,0])
    
    sr2 = 0.6/np.sqrt(2)
    x,y,z = sr2,sr2,0

    x = 0.6*np.cos(theta)
    y = 0.6*np.sin(theta)

    xcon,ycon = axisLength*x/0.6 , axisLength*y/0.6

    arc_line = np.transpose( [ [0,0,0], [x,y,0] ] )
    arc_line_cont = np.transpose( [[x,y,0] , [xcon,ycon,0]] )

    zt = np.sqrt(3)/2
    ver_line = np.transpose( [ [x,y,0], [x,y,zt] ] )
    ax.plot(ver_line[0], ver_line[1], ver_line[2], c=auxColor, linestyle=':')
    ax.scatter(x,y,zt,c=auxColor,s=40)

    ax.plot(arc_line[0], arc_line[1], arc_line[2], c='red', linestyle='-')
    ax.plot(arc_line_cont[0], arc_line_cont[1], arc_line_cont[2], c='grey', linestyle='-')
    ax.plot(arcx, arcy, arcz, c='blue', linewidth=1, linestyle='-')

    ax.scatter(x,y,z,c='tab:red',s=40)
    ax.add_collection3d(surface)

    ax.text(1.2,0.6,0,r'$\theta$', fontweight='bold', fontsize='large', color = 'blue', horizontalalignment='center', verticalalignment='center')
    ax.text(.7,.83,1,'(z)', fontsize='large', color = auxColor, horizontalalignment='center', verticalalignment='center')
    ax.text(.45,.2,0,'r', fontsize='large', color = 'red', horizontalalignment='center', verticalalignment='center')

    return

def plate(ax) :
    surface = s3d.PlanarSurface(4,facecolor=fcolor0)
    surface.set_edgecolor([0,0,0,0])
    
    sr2 = 0.75/np.sqrt(2)
    x,y,z = sr2,sr2,0

    x_line = np.transpose([ [0,y,0], [x,y,0] ])
    y_line = np.transpose([ [x,0,0], [x,y,0] ])
    ax.plot(x_line[0], x_line[1], x_line[2], c='red', linestyle='-')
    ax.plot(y_line[0], y_line[1], y_line[2], c='blue', linestyle='-')

    zt = np.sqrt(3)/2
    ver_line = np.transpose( [ [x,y,0], [x,y,zt] ] )
    ax.plot(ver_line[0], ver_line[1], ver_line[2], c=auxColor, linestyle=':')
    ax.scatter(x,y,zt,c=auxColor,s=40)

    ax.scatter(x,y,z,c='tab:red',s=40)
    ax.add_collection3d(surface)

    ax.text(.8,0.3,0,'y', fontsize='large', color = 'blue', horizontalalignment='center', verticalalignment='center')
    ax.text(.7,.83,1,'(z)', fontsize='large', color = auxColor, horizontalalignment='center', verticalalignment='center')
    ax.text(.3,.7,0,'x', fontsize='large', color = 'red', horizontalalignment='center', verticalalignment='center')

    return

# Construct figure, add surface, and plot ......................
funcArr = [plate, disk, cylinder, sphere]
minmax, ticks = (-1,1) , [-1,0,1]

for i in range(4) :
    fig = plt.figure(figsize=plt.figaspect(1))
    ax = plt.axes(projection='3d')
    ax.set(xlim=minmax, ylim=minmax, zlim=minmax)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_zticks(ticks)
    s3d.setupAxis(ax,offset=0.0,negaxis=False)
    ax.view_init(20, 25)
    funcArr[i](ax)

plt.show()