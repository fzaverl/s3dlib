import numpy as np
from matplotlib import pyplot as plt
import s3dlib.surface as s3d

#.. Guides: Vector Fields - Vectors at a Coordinate
#   uvw coordinates for cylindrical and spherical coordinates.

arc_seg = 25
theta = np.pi/3
phi = np.pi/3.5

thetaArc = np.linspace(0.0,theta, arc_seg)
arcz = np.linspace(0, 0, arc_seg)
r = arcz + 1.0
arcx = r * np.cos(thetaArc)
arcy = r * np.sin(thetaArc)

axisLength = 1.5
fcolor = [.7,.7,.6,0.2]
shade = .3

def uvw(ax,zoff=0) :
    """Draw uvw and xyz coordinate axes."""
    phiuvx = phi
    if zoff is not 0: phiuvx=np.pi/2
    z = np.cos(phiuvx)
    x = np.sin(phiuvx)*np.cos(theta)
    y = np.sin(phiuvx)*np.sin(theta)

    uvw_length = 0.5

    v_ux = uvw_length*x
    v_uy = uvw_length*y
    v_uz = uvw_length*z

    v_vx = -uvw_length*np.sin(theta)
    v_vy =  uvw_length*np.cos(theta)
    v_vz = 0.0 

    v_wx = -uvw_length*np.cos(phiuvx)*np.cos(theta)
    v_wy = -uvw_length*np.cos(phiuvx)*np.sin(theta)
    v_wz = uvw_length*np.sin(phiuvx)

    uvw_color = 'black'
    ax.quiver(x, y, z+zoff, v_ux  , v_uy  , v_uz,   color=uvw_color)
    ax.quiver(x, y, z+zoff, v_vx  , v_vy  , v_vz,   color=uvw_color)
    ax.quiver(x, y, z+zoff, v_wx  , v_wy  , v_wz,   color=uvw_color)
    
    tx_off = 1.2

    tu_x = x + tx_off*v_ux
    tu_y = y + tx_off*v_uy
    tu_z = z + tx_off*v_uz + zoff

    tv_x = x + tx_off*v_vx
    tv_y = y + tx_off*v_vy
    tv_z = z + tx_off*v_vz + zoff

    tw_x = x + tx_off*v_wx
    tw_y = y + tx_off*v_wy
    tw_z = z + tx_off*v_wz + zoff

    ax.text(tu_x,tu_y,tu_z,'u', fontsize='large', color = uvw_color,
            horizontalalignment='center', verticalalignment='center')
    ax.text(tv_x,tv_y,tv_z,'v', fontsize='large', color = uvw_color,
            horizontalalignment='center', verticalalignment='center')
    ax.text(tw_x,tw_y,tw_z,'w', fontsize='large', color = uvw_color,
            horizontalalignment='center', verticalalignment='center')

    #.. plot xyz axis origin...
    ax.plot([1,0,0], [0,0,1], [0,0,0], c='grey', linewidth=2, linestyle='-')
    ax.plot([0,0], [0,0], [0,1], c='grey', linewidth=2, linestyle='-')
    s3d.setupAxis(ax,offset=1.0,negaxis=False)
    ax.view_init(20, 25)
    ax.set_axis_off()

    return

def sphere(ax) :
    """Draw construction lines and surface."""
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

def cylinder(ax,zoff) :
    """Draw construction lines and surface."""
    surface = s3d.CylindricalSurface(4,facecolor=fcolor)
    surface.shade(depth=shade, direction=[-1,1,1])
    surface.set_edgecolor([0,0,0,0])
    
    x = np.cos(theta)
    y = np.sin(theta)
    z = zoff

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

# Construct figure, add surface, and plot ......................

fig = plt.figure(figsize=(7,3.5))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')
ax1.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1))
ax2.set(xlim=(-1,1), ylim=(-1,1), zlim=(-1,1))

cylinder(ax1,0.7)
uvw(ax1,0.7)
sphere(ax2)
uvw(ax2)

fig.tight_layout()
plt.show()