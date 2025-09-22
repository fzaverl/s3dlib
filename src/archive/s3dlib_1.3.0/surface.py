# Copyright (C) Frank Zaverl, Jr.
# See file LICENSE for license information.

from s3dlib import __version__

import math
import warnings
import copy
from functools import reduce
#from time import time

import numpy as np
from scipy import interpolate, spatial
from skimage import measure

import matplotlib as mpl
from matplotlib import cm, colors, image, tri, rcParams
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

#.. simple warning string
warnings.formatwarning = lambda m,c,n,l,line=0 : str(m)+'\n'

_MAXREZ = 8          # in subclass surfaces, maximum recursive triangulation.
_MAXPLREZ = 10       # in ParametricLine lines, maximum recursive segmentation.
_MINRAD = 0.01       # in subclassed polar and spherical surface split basetypes.
_SPLITSIZE = 0.0001  # for split grid geometries, phi offset from 2*pi
_MAXIEFACE = 250     # max of multi-type polyhedron faces to calc initedges
                     # ( to avoid non-optimize edge calc on large input messhes )

_COOR_KWARGS = {     # used for clipping & contour lines and line-filled surfaces
    'coor': [ [0], [0, 'planar','xyz','p','P'],
                   [1, 'c', 'C', 'cylinder',  'cylindrical'],
                   [2, 's', 'S', 'sphere', 'spherical'],
                   [3, 'x', 'X'],
                   [4, 'y', 'Y'],
                   [5, 'z', 'Z'] ]
}

_COORSYS = { "XYZ": 0, "POLAR": 1, "CYLINDRICAL": 2, "SPHERICAL": 3, "LONGLAT": 4 }

_ILLUM = (1,0,1)  # default light direction, (cmap_normals, shade, hilite & line direction).
_DFT_VLIM =     [ -1.0, 1.0 ]       # default face color value bounds
_DFT_VERTVLIM = [  0.0, 1.0 ]       # default vector magnitude value bounds
_DFT_ALR = 0.25   # default axis length ratio, head size to vector magnitude.
_DFT_VIEW = [  30, -60 ]            # default, axes elev and azim
_DFT_CMAP = rcParams['image.cmap']  # default, colormap

# +----------------------------------------------------------+
# |  S3Dlib v_1.3.0 was developed using Matplotlib v_3.8.0   |
# +----------------------------------------------------------+

# FutDev: Future Development notes.
# MROI:   Minimum Return On coding effort Investment.
#         (the juice isn't worth the squeeze for current release)

class Surface3DCollection(Poly3DCollection):
    """
    Base class for 3D surfaces.
    """
    _geomType = _COORSYS['XYZ']

    @staticmethod
    def chull(points, **kargs) :
        """
        Surface3DCollection Convex Hull for a set of 3D points. 
        
        Parameters
        ----------
        points : N x 3 float array-like

        Returns
        -------
        Surface3DCollection object

        """
        points = np.array(points)
        hull = spatial.ConvexHull(points)
        verts = points[hull.vertices]
        cPnt = verts.mean(axis=0)
        # set the faceIndices referencing the verts, not the points....
        revlkup = { val:i for i,val in enumerate(hull.vertices) }
        fIndices = [None]*len( hull.simplices )
        for i,indx in enumerate(hull.simplices):
            a,b,c = indx
            fIndices[i] = [ revlkup[a] , revlkup[b] , revlkup[c] ] 
        # set the face indices using the RHR for outward normals from the center.
        temp = [None]*len(fIndices)
        for i,findx in enumerate(fIndices) :
            a,b,c = findx
            fcenter = (verts[a] + verts[b] + verts[c])/3
            nordir = np.cross(verts[b]-verts[a],verts[c]-verts[b])
            dotprod = np.dot(fcenter-cPnt,nordir)
            temp[i] = [a,b,c] if dotprod>0 else [c,b,a]
        fIndices = temp
        surface = Surface3DCollection(verts,fIndices, **kargs)
        return surface

    @staticmethod
    def _get_cloud_points(N,dm) :
        """
        Grid of xyz mesh coordinates within a domain.

        Note: Cloud points only calculated once for surface sets in the same domain.
        """
        XD,YD,ZD = dm[0],dm[1],dm[2]
        dmT = dm.T
        A = dmT[1]-dmT[0]
        B=np.array( [ XD[0],YD[0],ZD[0] ])
        Nj = (N+1)*1j
        xyz = np.mgrid[XD[0]:XD[1]:Nj,YD[0]:YD[1]:Nj, ZD[0]:ZD[1]:Nj ]
        return xyz,A,B

    @staticmethod
    def _get_VF(cloud,sval,A,B) :
        """
        Surface vertices/face indices for an f(x,y,z) = constant.

        The cloud contains the scalar values.
        The surface of constant value, sval.
        A,B are specific to the cloud domain.
        """
        M = cloud.shape[0] - 1
        N = cloud.shape[1] - 1
        P = cloud.shape[2] - 1
        spacing = (1/M, 1/N, 1/P)
        level = 0.001    # FutDev: level should not be a constant.
        vol = cloud - sval
        verts, fvIndices, _, _ = measure.marching_cubes(vol, level, spacing=spacing)
        C = np.empty([len(verts),3])
        C[:] = A
        verts = C*verts + B
        return verts,fvIndices

    @staticmethod
    def _impl_cloud_surf(operation, drez, cloud ,domain, fval, name, **kwargs) :

        dm = _interpret_domain(domain)   # update for v_1.3.0
        XD,YD,ZD = dm[0],dm[1],dm[2]
        dmT = dm.T
        A = dmT[1]-dmT[0]
        B=np.array( [ XD[0],YD[0],ZD[0] ])
        if cloud is None :
            N = int( 10*drez )
            Nj = (N+1)*1j
            xyz = np.mgrid[XD[0]:XD[1]:Nj,YD[0]:YD[1]:Nj, ZD[0]:ZD[1]:Nj ]
            cloud = operation(xyz)

        verts, faces = Surface3DCollection._get_VF(cloud,fval,A,B)
        surface = Surface3DCollection( verts, faces, **kwargs)
        if name is not None: surface._geomName = name
        return surface

    @staticmethod
    def _impl_cloud_surfset(operation, drez, cloud ,domain, numb, name, **kwargs) :
    
        dm = _interpret_domain(domain)   # update for v_1.3.0
        XD,YD,ZD = dm[0],dm[1],dm[2]
        dmT = dm.T
        A = dmT[1]-dmT[0]
        B=np.array( [ XD[0],YD[0],ZD[0] ])
        if cloud is None :
            N = int( 10*drez )
            Nj = (N+1)*1j
            xyz = np.mgrid[XD[0]:XD[1]:Nj,YD[0]:YD[1]:Nj, ZD[0]:ZD[1]:Nj ]
            cloud = operation(xyz)
        # extract cmap or color from kwargs.
        kargDft = { 'color':None, 'cmap':_DFT_CMAP  }
        KW = KWprocessor(kargDft,'implsurfSet',**kwargs)
        userColor = KW.getVal('color')
        useDefColor = userColor is not None   # use color arg  if color is defined.
        cmapobj = KW.getVal('cmap')
        kwargs = KW.filter()
        if isinstance(cmapobj,str) : 
            cmapobj = mpl.colormaps[cmapobj]      # update for v_1.2.0

        delta = 1/(numb + 1)
        minVal, maxVal = np.min(cloud), np.max(cloud)
        rngVal = maxVal - minVal

        surface = None
        for i in np.linspace(delta, 1, numb, endpoint=False) :
            fval = minVal + i*rngVal
            verts, faces = Surface3DCollection._get_VF(cloud,fval,A,B)
            color = userColor if useDefColor else cmapobj(i)
            surf = Surface3DCollection( verts, faces, color=color, **kwargs)
            surface = surf if surface is None else surface + surf
        
        # following for a colorbar view
        surface.set_cmap(cmapobj)
        # note: min/max is cloud values, not surface range.
        surface.bounds['vlim'] = [ minVal, maxVal ]
        if name is not None: surface._geomName = name    # note: for this special case:
        if name is not None: surface.cname = name        #       cname = name

        return surface

    @staticmethod
    def implsurf(operation, drez=2.0, domain=1, fval=0.0, name=None, **kwargs) :
        """
        Surface3DCollection defined by an implicit function f(x,y,z).

        Parameters
        ----------
        operation : function object
            Implicit function that takes one
            argument, a 3xN Numpy array of x,y,z coordinates.
            The function returns a single scalar value.
            The value of zero defines the implicit surface.

        drez : Float, optional, default: 2.0
            Multiplier for the number of subdivisions within the domain.
            Frez values range from 1 to 10.  Function evaluation will be
            the order of frez cubed.

        domain : a number, list or array, default: 1
            The domain of the function evaluation.  For a number, n,
            the x,y,z axes domains will be [-n,n]. For a 1-dimensional 2-element list,
            [a,b] will be assigned the domain for all 3 axes.  Using a list of list (array),
            as [ [a,b],[c,d],[e,f] ], the domain is assigned individually for each of the
            three coordinate axes.

        fval : a number, default: 0.0
            Scalar value of the function at the surface.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Returns
        -------
        Surface3DCollection object

        """

        if not isinstance(drez, (int,float)) : 
            raise ValueError("Incorrect drez argument datatype passed to impf:", drez)
        if drez<1 or drez>10 :
            raise ValueError("Error: incorrect value drez passed to impf:", drez)

        surface = Surface3DCollection._impl_cloud_surf \
            (operation, drez,None,domain,fval,name,**kwargs)
        if name is None: name = _getFunctionName(operation,name)
        if name is not None: surface._geomName = name
        return surface

    @staticmethod
    def cloudsurf(cloud, domain=1, fval=0.0, name=None, **kwargs) :
        """
        Surface3DCollection defined for a constant value in a point cloud.

        Parameters
        ----------
        cloud : a (M,N,P) shaped array of point cloud values.

        domain : a number, list or array, default: 1
            The domain of the function evaluation.  For a number, n,
            the x,y,z axes domains will be [-n,n]. For a 1-dimensional 2-element list,
            [a,b] will be assigned the domain for all 3 axes.  Using a list of list (array),
            as [ [a,b],[c,d],[e,f] ], the domain is assigned individually for each of the
            three coordinate axes.

        fval : a number, default: 0.0
            Scalar value of surface contour within the cloude.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Returns
        -------
        Surface3DCollection object

        """

        surface = Surface3DCollection._impl_cloud_surf \
            (None,None,cloud,domain,fval,name,**kwargs)
        return surface

    @staticmethod
    def implsurfSet(operation, drez=2.0, domain=1, numb=2, name=None, **kwargs) :
        """
        Set of surface contours implicit function f(x,y,z) within a domain

        Parameters
        ----------
        operation : function object
            Implicit function that takes one
            argument, a 3xN Numpy array of x,y,z coordinates.
            The function returns a single scalar value.
            The value of zero defines the implicit surface.

        drez : Float, optional, default: 2.0
            Multiplier for the number of subdivisions within the domain.
            Frez values range from 1 to 10.  Function evaluation will be
            the order of frez cubed.

        domain : a number, list or array, default: 1
            The domain of the function evaluation.  For a number, n,
            the x,y,z axes domains will be [-n,n]. For a 1-dimensional 2-element list,
            [a,b] will be assigned the domain for all 3 axes.  Using a list of list (array),
            as [ [a,b],[c,d],[e,f] ], the domain is assigned individually for each of the
            three coordinate axes.

        numb : an integer, default: 2
            Number of surface contours within the domain.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Returns
        -------
        Surface3DCollection object

        """

        if not isinstance(drez, (int,float)) : 
            raise ValueError("Incorrect drez argument datatype passed to impf:", drez)
        if drez<1 or drez>10 :
            raise ValueError("Error: incorrect value drez passed to impf:", drez)

        surface = Surface3DCollection._impl_cloud_surfset \
            (operation, drez,None,domain,numb,name,**kwargs)
        if name is None: name = _getFunctionName(operation,name)
        if name is not None: surface._geomName = name
        return surface

    @staticmethod
    def cloudsurfSet(cloud, domain=1, numb=2, name=None, **kwargs) :
        """
        Surface3DCollection defined for a constant value in a point cloud.

        Parameters
        ----------
        cloud : a (M,N,P) shaped array of point cloud values.

        domain : a number, list or array, default: 1
            The domain of the function evaluation.  For a number, n,
            the x,y,z axes domains will be [-n,n]. For a 1-dimensional 2-element list,
            [a,b] will be assigned the domain for all 3 axes.  Using a list of list (array),
            as [ [a,b],[c,d],[e,f] ], the domain is assigned individually for each of the
            three coordinate axes.

        numb : an integer, default: 2
            Number of surface contours within the domain.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Returns
        -------
        Surface3DCollection object

        """

        surface = Surface3DCollection._impl_cloud_surfset \
            (None,None,cloud,domain,numb,name,**kwargs)
        return surface

    @staticmethod
    def _init_triangulateBase(rez,baseVcoor,baseFaceVertexIndices, midVectFunc) :
        """
        Recursively subdivide triangles.

        Each recursion subdivides each triangle by four.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.

        baseVcoor : V x 3 float array
            An array of 'v' number of xyz vertex coordinates.

        baseFaceVertextIndices : F x 3 int array
            An array of 'F' number of face vertex indices.

        midVectFunc : function object
            A function that takes two xyz coordinate (list of 3)
            representing a surface face edge.  Returns one xyz coordinate
            at the bisection of the edge, mapped on the surface.

        Returns
        -------
        indexObj : a dictionary of vertex indices, for a surface
            of F faces, E edges and V vertices.
            'face' : F x 3 int array
            'edge' : E x 2 int array

        vertCoor : V x 3 float array
            An array of V number of xyz vertex coordinates.

        """
        # triangulateBase to _triangulateBase,  update for v_1.3.0
        # _triangulateBase is now a instance method.
        # _init_triangulateBase only used for init of derived classes.
        # FutDev : use single method.

        vCoor = baseVcoor.copy()
        fvIndices = []
        evIndices = []
        #.....................................    
        def getVerticesLine(degree, leftIndex, rightIndex) :
            def recurs_edgeCoor(degree, indexOrigin, isRight, leftIndex, rightIndex, Eindex) :
                m = degree - 1
                delta = int(np.power(2,m))
                if isRight : delta = -delta
                i = indexOrigin - delta
                mid_coor = midVectFunc( vCoor[leftIndex], vCoor[rightIndex] )
                currentCoorIndex = len(vCoor)
                vCoor.append(mid_coor)
                Eindex[i] = currentCoorIndex
                if m <= 0 :
                    evIndices.append( [leftIndex, currentCoorIndex] )                    
                    evIndices.append( [currentCoorIndex, rightIndex] )                    
                    return
                recurs_edgeCoor(m,i,False,leftIndex,currentCoorIndex, Eindex)
                recurs_edgeCoor(m,i,True,currentCoorIndex,rightIndex, Eindex)
                return
            #.....................................    
            iOrigin = int(np.power(2,degree))
            Eindex = [None]*(iOrigin+1)
            Eindex[0] = leftIndex
            Eindex[iOrigin] = rightIndex
            if degree==0 :
                evIndices.append( [leftIndex, rightIndex] )
                return Eindex
            recurs_edgeCoor(degree, iOrigin, False, leftIndex, rightIndex, Eindex)
            return Eindex
        #.....................................    
        def recurs_faceIndices(A,B,C,rez,atCenter=True):
            if rez==0 :
                abc = [ A[0], B[0], C[0] ]               
                fvIndices.append( abc )
                if atCenter : evIndices.append( [ B[0], C[0] ] )
                return
            # -----------------------------------------------------
            J = int( (len(A)-1)/2 )
            K = J + 1
            # construct the center sub triangle edge vertices....
            if rez==1 :
                xz = [ A[1], C[1] ]
                yx = [ B[1], A[1] ]
                zy = [ C[1], B[1] ]
            else:               
                xz = getVerticesLine(rez-1, A[J], C[J] )
                yx = getVerticesLine(rez-1, B[J], A[J] )
                zy = getVerticesLine(rez-1, C[J], B[J] )
            # construct the 4 sub triangles...
            recurs_faceIndices(A[:K],xz,C[J:],rez-1)
            recurs_faceIndices(B[:K],yx,A[J:],rez-1)
            recurs_faceIndices(C[:K],zy,B[J:],rez-1)
            recurs_faceIndices(yx[::-1],zy[::-1],xz[::-1],rez-1,False)
            return
        #.....................................
        # generate the base edge indices for each face 
        baseFaceEdgeIndices = []
        edgeArrayIndexMatrix = [[None for x in range(len(vCoor))] for y in range(len(vCoor))]   
        for face in baseFaceVertexIndices :
            A = [ face[0],face[1]]
            B = [ face[1],face[2]]
            C = [ face[2],face[0]]
            edgeIndices = [A,B,C]
            baseFaceEdgeIndices.append(edgeIndices)
            for vertexIndex in edgeIndices :
                i = vertexIndex[0]
                j = vertexIndex[1]
                if edgeArrayIndexMatrix[i][j] is None:
                    edgeArrayIndexMatrix[i][j] = getVerticesLine(rez, i, j )
                    edgeArrayIndexMatrix[j][i] = edgeArrayIndexMatrix[i][j].copy()[::-1]
        # trianglulate the base triangles
        for face in baseFaceEdgeIndices :
            A = edgeArrayIndexMatrix[face[0][0]][face[0][1]]
            B = edgeArrayIndexMatrix[face[1][0]][face[1][1]]
            C = edgeArrayIndexMatrix[face[2][0]][face[2][1]]
            x = (rez > 0)
            recurs_faceIndices(A,B,C,rez, x)
        vertexCoor = np.array(vCoor) 
        indexObj = { 'face': np.array(fvIndices) , 'edge': np.array(evIndices) }
        return indexObj ,vertexCoor

    @staticmethod
    def _rectangulateBase(rez,baseVcoor,baseFaceVertexIndices, midVectFunc) :
        """
        Recursively subdivide quadrilaterals.

        Each recursion subdivides each quadrilateral by four.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the rectangle base faces.
            Rez values range from 0 to 7.

        baseVcoor : V x 3 float array
            An array of 'v' number of xyz vertex coordinates.

        baseFaceVertextIndices : F x 4 int array
            An array of 'F' number of face vertex indices.

        midVectFunc : function object
            A function that takes two xyz coordinate (list of 4)
            representing a surface face edge.  Returns one xyz coordinate
            at the bisection of the edge, mapped on the surface.

        Returns
        -------
        indexObj : a dictionary of vertex indices, for a surface
            of F faces, E edges and V vertices.
            'face' : F x 4 int array
            'edge' : E x 2 int array

        vertCoor : V x 3 float array
            An array of V number of xyz vertex coordinates.

        """
        # rectangulateBase to _rectangulateBase,  update for v_1.3.0
        vCoor = baseVcoor.copy()
        fvIndices = []
        evIndices = []
        #.....................................    
        def getVerticesLine(degree, leftIndex, rightIndex) :
            def recurs_edgeCoor(degree, indexOrigin, isRight, leftIndex, rightIndex, Eindex) :
                m = degree - 1
                delta = int(np.power(2,m))
                if isRight : delta = -delta
                i = indexOrigin - delta
                mid_coor = midVectFunc( vCoor[leftIndex], vCoor[rightIndex] )
                currentCoorIndex = len(vCoor)
                vCoor.append(mid_coor)
                Eindex[i] = currentCoorIndex
                if m <= 0 :
                    evIndices.append( [leftIndex, currentCoorIndex] )                    
                    evIndices.append( [currentCoorIndex, rightIndex] )                    
                    return
                recurs_edgeCoor(m,i,False,leftIndex,currentCoorIndex, Eindex)
                recurs_edgeCoor(m,i,True,currentCoorIndex,rightIndex, Eindex)
                return
            #.....................................    
            iOrigin = int(np.power(2,degree))
            Eindex = [None]*(iOrigin+1)
            Eindex[0] = leftIndex
            Eindex[iOrigin] = rightIndex
            if degree==0 :
                evIndices.append( [leftIndex, rightIndex] )
                return Eindex
            recurs_edgeCoor(degree, iOrigin, False, leftIndex, rightIndex, Eindex)
            return Eindex
        #.....................................    
        def recurs_faceIndices(A,B,C,D,rez):
            # similar to inner function of triangulateBase
            if rez==0 :
                abc = [ A[0], B[0], C[0], D[0] ]               
                fvIndices.append( abc )
                return
            # -----------------------------------------------------
            J = int( (len(A)-1)/2 )
            K = J + 1
            # midpoint indices ...............
            mid_coor = midVectFunc( vCoor[ D[J] ], vCoor[ B[J] ] )
            i_E = len(vCoor)
            vCoor.append(mid_coor)
            # arrays of the four inner quadrilaterals ..............
            pE = getVerticesLine(rez-1, A[J], i_E )
            qE = getVerticesLine(rez-1, B[J], i_E )
            rE = getVerticesLine(rez-1, C[J], i_E )
            sE = getVerticesLine(rez-1, D[J], i_E )
            # sudivide the four inner quadraterals .................
            recurs_faceIndices( A[:K], pE, sE[::-1], D[J:], rez-1)
            recurs_faceIndices( B[:K], qE, pE[::-1], A[J:], rez-1)
            recurs_faceIndices( C[:K], rE, qE[::-1], B[J:], rez-1)
            recurs_faceIndices( D[:K], sE, rE[::-1], C[J:], rez-1)
            return
        #.....................................

        # generate the base edge indices for each face 
        baseFaceEdgeIndices = []
        edgeArrayIndexMatrix = [[None for x in range(len(vCoor))] for y in range(len(vCoor))]   
        for face in baseFaceVertexIndices :
            A = [ face[0],face[1]]
            B = [ face[1],face[2]]
            C = [ face[2],face[3]]
            D = [ face[3],face[0]]
            edgeIndices = [A,B,C,D]
            baseFaceEdgeIndices.append(edgeIndices)
            for vertexIndex in edgeIndices :
                i = vertexIndex[0]
                j = vertexIndex[1]
                if edgeArrayIndexMatrix[i][j] is None:
                    edgeArrayIndexMatrix[i][j] = getVerticesLine(rez, i, j )
                    edgeArrayIndexMatrix[j][i] = edgeArrayIndexMatrix[i][j].copy()[::-1]
        # rectangulate the base quadrilaterals
        for face in baseFaceEdgeIndices :
            A = edgeArrayIndexMatrix[face[0][0]][face[0][1]]
            B = edgeArrayIndexMatrix[face[1][0]][face[1][1]]
            C = edgeArrayIndexMatrix[face[2][0]][face[2][1]]
            D = edgeArrayIndexMatrix[face[3][0]][face[3][1]]
            recurs_faceIndices(A,B,C,D,rez)
        vertexCoor = np.array(vCoor) 
        indexObj = { 'face': np.array(fvIndices) , 'edge': np.array(evIndices) }
        return indexObj ,vertexCoor

    @staticmethod
    def coor_convert(xyz, tocart=False) :
        """Coordinate transformation.
        
           To be overriddden by any subclass not in Cartesian coordinates.
        """
        xyz = np.array(xyz)
        return xyz

    @staticmethod
    def _map3Dto2Dsquare(xyz) :
        """ Map surface to rectangular unit coordinates (0,1). """

        x,y,_ = xyz
        a = (x+1)/2
        b = (y+1)/2
        return [a,b]

    @staticmethod
    def _sfev(rez,fevArray) :
        """ Calculate number of surface faces, edges, vertices. """
        # vestigial
        Fo, phi, Vx = fevArray
        eta = np.power(2,rez)
        F = eta*eta*Fo
        E = np.int( 3*F/2 + eta*phi )
        V = np.int( F/2 + eta*phi + Vx )
        return [F,E,V]

    @classmethod
    def _grid_face_indices(cls,U, V, isSplit=False, isTria=False) :
        """ method only called by the classmethod _square_net or
            the CylindricalSurface classmethod _cyl_vol_square_net
        """
        if isSplit :
            t = np.linspace(0,V*(U+1)-1,V*(U+1),dtype=int )
            ts = t.reshape( (V,U+1) )
            a = ts[:,0:U].flatten()
            findx = np.array( [ a, a+1, a+U+2, a+U+1 ] )
        else :    
            a = np.linspace(0,V*U-1,V*U,dtype=int )
            bs = a.reshape( (V,U) ).T
            order = np.append(np.linspace(1,U-1,U-1,dtype=int),0)
            b = bs[order].T.flatten()
            findx = np.array( [ a, b, b+U, a+U ] )
        if isTria :
            a,b,c,d = findx
            findx = np.concatenate( ( [a,b,c],[c,d,a] ), axis=1   )
        return findx.T

    @classmethod
    def _square_net(cls, numV, numU, basetype,minRad, splitSize, name, **kwargs) :
        """ method only called by inherited class classmethods 'grid'. """

        divLimits = [1,3]
        if cls._geomType is _COORSYS['XYZ']  : divLimits = [1,1]
        if cls._geomType is _COORSYS['SPHERICAL'] : divLimits = [2,3]
        if numV < divLimits[0] :
            raise ValueError('Grid division {} must be greater or equal to {}'.format(str(numV),divLimits[0] ))
        if numU < divLimits[1] :
            raise ValueError('Grid division {} must be greater or equal to {}'.format(str(numU),divLimits[1] ))

        if basetype == None : basetype = 'd'
        if minRad is None : minRad = _MINRAD
        if splitSize is None : splitSize = _SPLITSIZE

        basetype = basetype.lower()
        baselist = ['d','s','q','w','r','x']
        if basetype not in baselist :
            raise ValueError("Grid basetype '{}' is not recognized. Possible values are: {}.".format(basetype,list(baselist)))

        isSplit = (basetype == 's') or (basetype == 'w') or (basetype == 'x')
        needsTriag = not ( (basetype == 'r') or (basetype == 'x')  )
        hasCenter = (basetype == 'd') or (basetype == 's')
        # hasCenter does not apply for cylindrical grids of this type.
        if cls._geomType is _COORSYS['CYLINDRICAL'] : hasCenter=False

        # Note: numX and num_X are the number of divisions and the number of coor positions.
        if cls._geomType is _COORSYS['XYZ']  :
            isSplit = True
            hasCenter = False
            num_U,min_U, max_U = numU+1,-1.0,1.0
            num_V = numV + 1
            coor_U = np.linspace(min_U,max_U,num_U)
            U_coor = np.tile(coor_U,num_V)
        else :
            # T angle ranges, used for disk,cylinder and sphere.....
            num_U = numU+1 if isSplit else numU
            if hasCenter :
                num_V = numV if cls._geomType is _COORSYS['POLAR'] else numV-1
            else :
                num_V = numV+1
            min_U = 0.0
            max_U = (1-1/num_U)*2*np.pi
            if isSplit : max_U = (2-splitSize)*np.pi
            coor_U = np.linspace(min_U,max_U,num_U)
            U_coor = np.tile(coor_U,num_V)
        totVerts = num_U*num_V

        if cls._geomType is _COORSYS['XYZ'] :
            val_W = 0.0  # z is constant at 0.
            coor_V =  np.linspace(-1.0,1.0,num_V)
            order_abc = lambda u,v,w :[u,v,w]
        if cls._geomType is _COORSYS['POLAR'] :
            val_W = 0.0  # z is constant at 0.
            end_V = 1/numV if hasCenter else minRad
            coor_V = np.linspace(1.0,end_V,num_V)
            order_abc = lambda u,v,w : [v,u,w]
        if cls._geomType is _COORSYS['CYLINDRICAL'] :
            val_W = 1.0  # r is constant at 1.
            coor_V =  np.linspace(-1.0,1.0,num_V)
            order_abc = lambda u,v,w : [w,u,v]
        if cls._geomType is _COORSYS['SPHERICAL'] :
            val_W = 1.0  # r is constant at 1.
            st_angle = np.pi/numV if hasCenter else np.arcsin(minRad)
            stt_V = np.pi - st_angle
            end_V = st_angle
            coor_V =  np.linspace(stt_V,end_V,num_V)
            order_abc = lambda u,v,w : [w,u,v]

        W_coor = np.tile(val_W,totVerts)
        V_coor = np.tile(coor_V,(num_U,1))
        V_coor = np.ravel(V_coor,order='F')

        abc = order_abc(U_coor,V_coor,W_coor)
        abc = np.array(abc)
        xyz = cls.coor_convert(abc,True)
        verts = np.array(xyz).T
        numVt = num_V-1 if hasCenter else numV
        faceIndices = cls._grid_face_indices(numU,numVt,isSplit,needsTriag)

        if hasCenter :
            # add top triangular faces for polar & spherical surfaces ...
            topVert = np.array([[0.0,0.0,val_W]])
            topIndex = totVerts
            startIndex = topIndex - num_U
            endIndex = topIndex-2 if isSplit else topIndex-1
            a = np.linspace(startIndex,endIndex,numU,dtype=int)
            b = a+1 if isSplit else np.append( a[1:], [startIndex] )
            c = np.full(numU,topIndex)
            abc = np.array([a,b,c]).T
            verts = np.append(verts,topVert, axis=0)
            faceIndices = np.append(faceIndices,abc, axis=0)
            if cls._geomType is _COORSYS['SPHERICAL'] :
                # add bottom triangular faces for spherical surfaces ...
                btmVert = np.array([[0.0,0.0,-val_W]])
                btmIndex = totVerts+1
                startIndex,endIndex = 0,numU-1
                a = np.linspace(startIndex,endIndex,numU,dtype=int)
                b = a+1 if isSplit else np.append( a[1:], [startIndex] )
                c = np.full(numU,btmIndex)
                bac = np.array([b,a,c]).T  # bac instead for abc for outward normals.
                verts = np.append(verts,btmVert, axis=0)
                faceIndices = np.append(faceIndices,bac, axis=0)

        surface = cls(**kwargs)  # set to the class default object
        # reset verts, faceIndices and edges (retain class coordinate system).
        surface.vertexCoor = verts
        surface.fvIndices = faceIndices
        surface.evIndices = surface._get_edges_from_faces(faceIndices)
        surface.vfIndicesList = surface._vertexCommonFaceIndices()
        surface.efIndicesList = surface._edgeCommonFaceIndices()
        surface.set_verts( verts[faceIndices] )
        surface._set_geometric_bounds()
        surface.name = name
        surface._basetype = 'grid_'+basetype
        return surface

    @classmethod
    def _fev(cls,rez,basetype,dercls) :
        """ Sets up _sfev method, based on calls from the subclass dercls. """
        # vestigial
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None: basetype = dercls._default_base
        basesurf = dercls._base_surfaces
        if basetype not in basesurf :
            raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(basesurf)))
        return cls._sfev(rez, basesurf[basetype]['fevParam'] )

    @classmethod
    def _uvw_to_xyz(cls,rst,abc) :
        """ Rotational transformation at a coordinate.
            To be overridden by any subclass not in Cartesian coordinates.
        """
        return np.transpose(rst)

    @classmethod
    def _disp_to_xyz( cls, uvw, theta=None, phi=None) :
        """ Rotational transformation from rotation angles. """

        if theta is None : return uvw
        # each set of uvw has a different rotation matrix
        # that must be evaluated one at a time. MROI
        if phi is None : phi = np.zeros( len(theta) )
        uvwT = np.transpose(uvw)
        dispVect = []
        for i in range(len(theta)) :
            rotMx = eulerRot(theta[i],phi[i], inrad=True)
            vec = np.dot( uvwT[i], rotMx )
            dispVect.append(vec)
        return dispVect

    def _triangulateBase(self,rez,midVectFunc) :
        """
        Recursively subdivide triangles.

        Each recursion subdivides each triangle by four.

        Includes updated calculations for vertex normals and,
        if vertex values are not None, vertex values.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.

        midVectFunc : function object
            A function that takes two xyz coordinate (list of 3)
            representing a surface face edge.  Returns one xyz coordinate
            at the bisection of the edge, mapped on the surface.

        Returns
        -------
        indexObj : a dictionary of vertex indices, for a surface
            of F faces, E edges and V vertices.
            'face' : F x 3 int array
            'edge' : E x 2 int array

        vertCoor : V x 3 float array
            An array of V number of xyz vertex coordinates.

        vertexValues V float array
            An array of values for each vertex.
        """
        # update 1.2.0 to 1.3.0 :
        # this object method replaces the class method triangulateBase.

        # NOTE: arrays are converted to lists (mutable) due to scope 
        #       accessability within the nested inner functions.
        
        # initialize vertex values, coordinates, phong_normals. Also
        #            face vertex-indicies and edge vertex-indicies.
        vCoor = list(self.vertexCoor)
        pvNrm = list(self._get_vertex_normals().copy())
        vVals = None
        if self.vertexVals is not None :
            vVals = list(self.vertexVals)
        fvIndices = []
        evIndices = []

        baseFaceVertexIndices = self.fvIndices
        #.....................................    
        def getVerticesLine(degree, leftIndex, rightIndex) :
            def recurs_edgeCoor(degree, indexOrigin, isRight, leftIndex, rightIndex, Eindex) :
                m = degree - 1
                delta = int(np.power(2,m))
                if isRight : delta = -delta
                i = indexOrigin - delta
                mid_coor = midVectFunc( vCoor[leftIndex], vCoor[rightIndex] )
                currentCoorIndex = len(vCoor)
                vCoor.append(mid_coor)
                # .... update for v_1.3.0
                if vVals is not None :
                    mid_val = (vVals[leftIndex] + vVals[rightIndex])/2
                    vVals.append(mid_val)
                
                mid_norm = pvNrm[leftIndex] + pvNrm[rightIndex]
                mid_norm = np.divide( mid_norm, np.linalg.norm(mid_norm) )
                pvNrm.append(mid_norm)
                # .......................
                Eindex[i] = currentCoorIndex
                if m <= 0 :
                    evIndices.append( [leftIndex, currentCoorIndex] )                    
                    evIndices.append( [currentCoorIndex, rightIndex] )                    
                    return
                recurs_edgeCoor(m,i,False,leftIndex,currentCoorIndex, Eindex)
                recurs_edgeCoor(m,i,True,currentCoorIndex,rightIndex, Eindex)
                return
            #.....................................    
            iOrigin = int(np.power(2,degree))
            Eindex = [None]*(iOrigin+1)
            Eindex[0] = leftIndex
            Eindex[iOrigin] = rightIndex
            if degree==0 :
                evIndices.append( [leftIndex, rightIndex] )
                return Eindex
            recurs_edgeCoor(degree, iOrigin, False, leftIndex, rightIndex, Eindex)
            return Eindex
        #.....................................    
        def recurs_faceIndices(A,B,C,rez,atCenter=True):
            if rez==0 :
                abc = [ A[0], B[0], C[0] ]               
                fvIndices.append( abc )
                if atCenter : evIndices.append( [ B[0], C[0] ] )
                return
            # -----------------------------------------------------
            J = int( (len(A)-1)/2 )
            K = J + 1
            # construct the center sub triangle edge vertices....
            if rez==1 :
                xz = [ A[1], C[1] ]
                yx = [ B[1], A[1] ]
                zy = [ C[1], B[1] ]
            else:               
                xz = getVerticesLine(rez-1, A[J], C[J] )
                yx = getVerticesLine(rez-1, B[J], A[J] )
                zy = getVerticesLine(rez-1, C[J], B[J] )
            # construct the 4 sub triangles...
            recurs_faceIndices(A[:K],xz,C[J:],rez-1)
            recurs_faceIndices(B[:K],yx,A[J:],rez-1)
            recurs_faceIndices(C[:K],zy,B[J:],rez-1)
            recurs_faceIndices(yx[::-1],zy[::-1],xz[::-1],rez-1,False)
            return
        #.....................................
        # generate the base edge indices for each face 
        baseFaceEdgeIndices = []
        edgeArrayIndexMatrix = [[None for x in range(len(vCoor))] for y in range(len(vCoor))]   
        for face in baseFaceVertexIndices :
            A = [ face[0],face[1]]
            B = [ face[1],face[2]]
            C = [ face[2],face[0]]
            edgeIndices = [A,B,C]
            baseFaceEdgeIndices.append(edgeIndices)
            for vertexIndex in edgeIndices :
                i = vertexIndex[0]
                j = vertexIndex[1]
                if edgeArrayIndexMatrix[i][j] is None:
                    edgeArrayIndexMatrix[i][j] = getVerticesLine(rez, i, j )
                    edgeArrayIndexMatrix[j][i] = edgeArrayIndexMatrix[i][j].copy()[::-1]
        # trianglulate the base triangles
        for face in baseFaceEdgeIndices :
            A = edgeArrayIndexMatrix[face[0][0]][face[0][1]]
            B = edgeArrayIndexMatrix[face[1][0]][face[1][1]]
            C = edgeArrayIndexMatrix[face[2][0]][face[2][1]]
            x = (rez > 0)
            recurs_faceIndices(A,B,C,rez, x)
        vertexCoor = np.array(vCoor)
        vertexValues = np.array(vVals) 
        phongNorms = np.array(pvNrm)
        #print('[_triangulateBase]  shape of vertexCoor,phongNorms',vertexCoor.shape,phongNorms.shape)
        indexObj = { 'face': np.array(fvIndices) , 'edge': np.array(evIndices) }
        return indexObj ,vertexCoor, vertexValues, phongNorms

    def _postProc_surfaceColors(self,onlyFaces=False) :
        """ Matplotlib will 'fillin' colors during rendering using the
            facecolor array. (this is NOT the assigned color to faces).
            However, S3Dlib will use the inherited _facecolors 
            for the assigned face colors.
            Whenever color is 'externally' assigned or reassigned,
            need to directly assign colors for each face for any
            method that uses facecolors.  This occurs initially
            in the method before colors are assigned            
            (set_color, +, clipping, shading, edges, etc.)
            
            Note that when initialization using **kargs,
            onlyFaces should be True so don't reset edgecolors
            if they are assigned.
        """
        N = len(self.fvIndices)      # number of faces
        initFC = self._dfacecolors   # update for v_1.2.0
        inLen = len(initFC)          # number of facecolors
        if inLen == N : return
        if isinstance(initFC,np.ndarray) : initFC = initFC.tolist()
        n = math.ceil(N/inLen)
        colorcoll = initFC*n
        if onlyFaces : self.set_facecolor(colorcoll[0:N])
        #if onlyFaces : self._dfacecolors = colorcoll[0:N]   # update for v_1.2.0
        else :         self.set_color(colorcoll[0:N])
        return

    def _set_index_relations(self) :
        """ from the vertex and faceIndices, set the indices:
            evIndices:     Array, for each edge, 2 vertex indicies.
            vfIndicesList: List,  for each vertex, face indices.
            efIndicesList: List,  for each edge, 0 to 2 face indices.
        """
        edgeIndices = self._get_edges_from_faces(self.fvIndices)          
        self.evIndices = np.array(edgeIndices)
        self.vfIndicesList = self._vertexCommonFaceIndices()
        self.efIndicesList = self._edgeCommonFaceIndices()
        return

    @staticmethod
    def _extract_vals_from_inputCoor(ponts) :
        """
        Separate (N,4) into (N,3) and (N) arrays,

        or just (N,3) and None.
        Used by __init__ and subclassed classmethods pntsurf
        """
        # update for v_1.3.0
        tempCoor = np.array(ponts,dtype=float)
        vVals = None
        if tempCoor.shape[1] == 3 :
            vCoor = tempCoor
        else :
            vCoor = tempCoor[:,:3]
            vVals = tempCoor[:,3]
        return vCoor,vVals

    def __init__(self, vertexCoor, faceIndices, name=None, vcolor=None, **kwargs) :
        """ 
        3D surface of connected F faces, each face with N number of vertices.
        The surface has V number of vertices and E number of edges.
        
        Parameters
        ----------
        vertexCoor : V x 3, or V X 4, float array-like 
            An array of V number of xyz vertex coordinates.
            If V X 4, scaler values are assigned to each vertex.

        faceIndices : F x N int list or array
            A list or array of F number of faces with N vertex indices.
            For lists, N may be different for each face.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        vcolor : V x 3 float array-like, default: None.
            An array of V number of vertex colors.
            (not implemented, reserved for future development)

        Other Parameters
        ----------------
        **kwargs
            All other arguments are passed on to mpl_toolkits.mplot3d.art3d.Poly3DCollection.
            Valid argument names include: color, edgecolors, facecolors, linewidths, cmap.

        """

        # following for a 'control limit' parameters in kwargs.
        KW = KWprocessor({'maxieface':_MAXIEFACE},'Surface3DCollection.__init__',False,**kwargs)
        maxieface = KW.getVal('maxieface')
        kwargs = KW.filter()

        dataPts = np.array(vertexCoor,dtype=float)
        # update for v_1.3.0
        dataPts, vVal = Surface3DCollection._extract_vals_from_inputCoor(dataPts)
        self.vertexCoor = dataPts
        self.vertexVals = vVal

        # for possible 'smooth' triangulated face shading.  # update for v_1.3.0
        self.vertexNorms = None        # geometric from face normals.
        self.phongNorms = None         # linear approximate on a flat face during triangulation.
        self._initedges = None         # assigned for irregulary type face meshes. (1.3 update)

        self.baseVertexCoor = np.array(self.vertexCoor,dtype=float)  # update for v_1.3.0
        self.fvIndices = np.array(faceIndices, dtype=object)  # update for v_1.2.0

        if self.fvIndices.ndim == 1 :
            # assume input is a list of faces with different number of edges.
            
            # the following restricts usage for non-optimized edge calc for large
            # meshes containing non-uniform polyhedral face types.  Otherwize, only
            # triangulated meshes are available.
            if  len(self.fvIndices) <= maxieface :    # update for v_1.3.0
                self._initedges = self._get_edges_from_faces(self.fvIndices)
                self._initverts = np.array(self.vertexCoor,dtype=float)  # update for v_1.3.0

            fvIndNew = self._triangulate_faceIndices(faceIndices)
            fcolor = None
            # note: if edgecolors are not defined, default to defined facecolor list,
            #       not the defined color ( continue with different 'magic' ) MROI.
            # following 'may' be vestigial to be removed ?
            if 'color' in kwargs : fcolor = kwargs['color']
            if 'facecolor' in kwargs : fcolor = kwargs['facecolor']
            if colors.is_color_like(fcolor) :
                        # single color assigned 
                fvIndNew = self._triangulate_faceIndices(faceIndices)
            else :
                if fcolor is None :
                        # no color assigned
                    fvIndNew = self._triangulate_faceIndices(faceIndices)
                else :  # list of colors assigned  
                    fvIndNew, faceColorsNew = self._triangulate_faceIndices_faceColors(faceIndices,fcolor)
                    kwargs['facecolor'] = faceColorsNew 
                    kwargs['color'] = faceColorsNew 
            self.fvIndices = np.array(fvIndNew)
        else :
            self.fvIndices = np.array(faceIndices)  # update for v_1.2.0
        self._set_index_relations()
        f = self.fvIndices
        v = self.vertexCoor
        super().__init__(v[f], **kwargs)
        # note: in following, only set the facecolors since Matplotlib auto-setting
        #       the edges, if not explicitely passed with edgecolor kwargs.
        self._postProc_surfaceColors(True) 

        #self._geomName = None    # Redundant from next line, internal tracking to determine if set...
        self.name = name
        self.valuesName = None
        # ========================================================================
        self.coorType = _COORSYS["XYZ"]
        self.baseSurfaceColor = np.array(self._dfacecolors[0]).flatten().tolist()  # update for v_1.2.0
       
        self.vertexColor = np.array(self._dfacecolors[0]).copy()[np.newaxis,:]     # update for v_1.2.0
        # following assigned when surface color applied using an operation
        self.vertexValues = None
        self.faceValues = None

        if len (self._dedgecolors) ==  0 :     # update for v_1.2.0
            self.baseEdgeColor = self.baseSurfaceColor 
        else :
            self.baseEdgeColor = np.array(self._dedgecolors[0]).flatten().tolist()     # update for v_1.2.0
        
        self.isSurfaceColored     = False       # flag set True for color array export.
        self.onlyColoredFaces     = False       # flag set True for cmap normals.
        self._surfaceShapeChanged = False       # control flag for image operations      
        self._isCompositeSurface  = False       # control flag for a composite object.

        self.disp_Vector = None                 # overridden in subclass.
        self._normal_scale = 1.0                # overridden in subclass default
        
        self.restoredClip = None                # dictionary to restore from clip operation.
        self.vector_field = None                # for vector fields when transformed.
        self._svd_dict    = None                # results from map_geom_from_svd
        self.scale = [1,1,1]                    # for axis when transformed.
        self.translation = [0,0,0]              # for axis when transformed.
        self.rotation =  np.identity(3)         # for axis when transformed.
        self._rez = 0                           # used by child classes.
        self._totRez = 0                        # account for triangulation, update for v_1.3.0

        self._triNumber           = 0           # flag, triangulation vertex normals
        self._bounds = {}
        self._set_geometric_bounds()
        self._bounds['vlim']     = _DFT_VLIM       
        self._bounds['vertvlim'] = _DFT_VERTVLIM 
    
    def __add__(self, othr) :
        """
        Combine two surface objects into a single surface object.

        Parameters
        ----------
        othr : Surface3DCollection object 

        Returns
        -------
        Surface3DCollection object

        """

        if not isinstance(othr, Surface3DCollection) :
            raise ValueError('Add operations can only apply between Surface3DCollection objects.')

        # MROI: should check if 2 surface faces are both triangular or quadrilateral

        topVerCoor = self.vertexCoor
        botVerCoor = othr.vertexCoor
        totVerCoor = np.append(topVerCoor,botVerCoor,axis=0)
        nextIndex = len(topVerCoor)

        topFVindices = self.fvIndices
        botFVindices = othr.fvIndices
        FVtail = np.add(botFVindices,nextIndex)
        totFVindices = np.append(topFVindices,FVtail,axis=0)

        obj = Surface3DCollection(totVerCoor,totFVindices)
        obj._isCompositeSurface = True
        
        top_onearr = np.ones( len(topFVindices) )[:, np.newaxis]
        top_orig_colors = self._dfacecolors  # update for v_1.2.0
        if len(top_orig_colors) == 1 : top_orig_colors = top_orig_colors*top_onearr

        bot_onearr = np.ones( len(botFVindices) )[:, np.newaxis]
        bot_orig_colors = othr._dfacecolors  # update for v_1.2.0
        if len(bot_orig_colors) == 1 : bot_orig_colors = bot_orig_colors*bot_onearr

        total_colors = np.append(top_orig_colors,bot_orig_colors,axis=0)
        obj.set_color(total_colors)        

        obj._set_geometric_bounds()
        return obj

    def __str__(self) :
        # parent class doesn't have specific properties...
        # Note: self._rez may be a string or int.  probably not good but MROI
        try:     bs = ' ( {} , {} )'.format(self._rez,self._basetype)
        except:  bs = ' '
        numFaces = len(self.fvIndices)
        numVerts = len(self.vertexCoor)
        name = self.__class__.__name__
        sz = ':  faces: {},  vertices: {}'.format(numFaces,numVerts)
        val = name+bs+sz
        return val

    def _calc_rez(self,N0) :
        # string val for equivalent grid rez 
        numFaces = len(self.fvIndices)
        cRez = np.log2(numFaces/N0)/2
        return '{:.1f}'.format(cRez)

    def _vertexCommonFaceIndices(self) :
        """ list of face indices at each vertex. """

        # note: a vertex may have NO faces (empty faceIndex list)
        #       due to clipping or a vertex not used for a face, so:
        vfIndicesList = [ [] for i in range(len(self.vertexCoor))]
        # note: usually vertices may have 1 to 6 common faces
        for i,face in enumerate(self.fvIndices) :
            for vInx in face :
                if vfIndicesList[vInx] is None :  vfIndicesList[vInx] = [i]
                else :                            vfIndicesList[vInx].append(i)
        return vfIndicesList
        
    def _edgeCommonFaceIndices(self) :
        """ list of face indices adjacent to each edge. """

        efIndicesList = [None]*len(self.evIndices)
        # note: edges may have 1 or 2 common faces
        for i,eIndex in enumerate(self.evIndices) :
            head_vfi = self.vfIndicesList[ eIndex[0] ]
            tail_vfi = self.vfIndicesList[ eIndex[1] ]
            for fi in head_vfi :
                if (fi in tail_vfi) :
                    if efIndicesList[i] is None:   efIndicesList[i] = [fi]  
                    else:                          efIndicesList[i].append(fi)
        return efIndicesList

    def _contourLineCollection( self, dotprod, projVect, *args ) :
        """
            ColorLine3DCollection of segments on faces interecting a plane.

            dotprod  - dot product between each vertex and plane normal.
            projVect - interection plane normal
            *args    - list of distances of the plane to the origin.
        """
        # FutDev: There are 'vestigal' loops during early development
        # debugging which need to be cleaned out. (MROI for release)
        # ---------------------------------------------------------------
        def getLineVals( dotprod,aLen,projVect) :
            # ................................................................
            def order_value( face,eIndices, conCoor, planeNormal) :
                A = eIndices[0]
                B = eIndices[1]
                fvi = self.fvIndices[face]
                verts = self.vertexCoor[ fvi ]
                faceNormal = self._get_face_normals([verts])

                edgeVert_A = conCoor[A]
                edgeVert_B = conCoor[B]
                A_to_B = np.subtract(edgeVert_B, edgeVert_A)
                edgeCross = np.cross(A_to_B, planeNormal )
                edgeCross= np.divide( edgeCross, np.linalg.norm(edgeCross) )
                test = np.dot(faceNormal, edgeCross)
                return [A,B] if test>0 else [B,A]
            # ................................................................
            
            # .. get intersect coor for each edge........
            evi = self.evIndices
            verts = self.vertexCoor
            dp_a = aLen - dotprod
            eElev = dp_a[evi]
            # eElev is the distance of the two edge vertices above/below the plane.
            # edges cross the plane if one vertex is above, the other below.
            # hence the product of the two distances is negative if edge crosses the plane.
            edChk = eElev[:,0]*eElev[:,1]
            xCoor = []  # intersect coordinates.
            evix = []   # the edge index of the intersect coordinate.
            for i in range( len(evi) ) :
                if edChk[i] < 0.0 :
                    p1 = verts[ evi[i,0] ]
                    p2 = verts[ evi[i,1] ]
                    dp1 = dotprod [ evi[i,0] ]
                    dp2 = dotprod [ evi[i,1] ]
                    a_mdp1 = dp_a [ evi[i,0] ]
                    lmb = a_mdp1/( dp2 - dp1 )
                    px = lmb*np.subtract(p2,p1) + p1
                    xCoor.append(px)
                    evix.append(i)

            # .. determine faces having edges which have an intersect.....
            if len(xCoor) == 0 : return None,None,None
            efi = self.efIndicesList
            numbXedges = len(evix)
            findexDict = {}  # faces(key=face index) which have two edges (value=xCoor index)
            for i in range(numbXedges) :
                ei = evix[i]
                efiArr = efi[ei]
                for fi in efiArr :
                    if fi in findexDict :
                        curVal = findexDict[fi]
                        if isinstance(curVal,list) :
                            # DevNote: this modifications is a consequence of quadrilateral
                            # faces, which may have more than one intersectioon to a plane.
                            # This only eliminates runtime error.  Solution is to use a surface
                            # with lrez > 0, then generated contours. 
                            newVal = [curVal[0],i]
                        else :
                            newVal = [curVal,i]                       
                        findexDict[fi] = newVal
                    else :
                        findexDict[fi] = i
            # remove face interset points if face only has one intersect edge
            # ( intersect may be on a vertex )
            #
            test = [ x for x in findexDict.values() if type(x) is list ]
            if len(test) == 0 : return None,None,None
            #
            face, coorInx = [], []
            for key,value in findexDict.items() :
                if type(value) is list :
                    valueX = order_value(key,value,xCoor,projVect)
                    coorInx.append(valueX)
                    face.append(key)

            if len(coorInx) == 0 : return None,None,None
            return np.array(xCoor), coorInx, face
        # ---------------------------------------------------------------
        def getColors(colors, indices) :
            con_colors = [None]*len(indices)
            for i in range(len(indices)) :
                con_colors[i] = colors[indices[i]]
            return con_colors
        # ---------------------------------------------------------------
        # get the surface colors which will be applied to the default contour colors
        self._postProc_surfaceColors()
        surf_colors = self._dfacecolors  # update for v_1.2.0

        if type(args[0]) is tuple : args=args[0]
        line = None
        for aLen in args :
            vertexCoor,segmIndices,faceIndices = getLineVals( dotprod,aLen,projVect  )
            if vertexCoor is None :
                warnings.warn( 'Surface contourLines not found for dist value of {}'.format(aLen) )
                continue
            tline = ColorLine3DCollection(vertexCoor,segmIndices)
            if tline is None : continue
            tline.set_color(getColors(surf_colors,faceIndices)) # default surface color
            if line is None: line = tline
            else: line += tline 
        return line

    def _transformVector(self, orgCoor, operation, returnxyz) :
        """ Functional transformation of coordinates. """
        xyz = np.transpose(orgCoor)
        abc = self.coor_convert(xyz)
        rst = np.array(operation(abc))
        if returnxyz : XYZ = rst
        else :         XYZ = self.coor_convert( rst , tocart=True )
        return np.transpose(XYZ)

    def _get_vectorfield_Line3DCollection(self, vector_field, color, width, alr) :
        """Vector3DCollection at surface vertices.""" 

        self.vector_field = vector_field
        
        lcol = Vector3DCollection(self.vertexCoor,vector_field,alr,colors=color, linewidths=width)
        # --- 1.1 mod --------------------
        lcol.coorType = self.coorType
        lcol.uvwOrientationType = self.coorType
        return lcol

    def _get_vertex_normals(self) :
        """Unit normals of at vertex coordinates. """

        # The normalized unweighted average of the surface
        # normals of the faces that contain that vertex.

        verts = self.vertexCoor[ self.fvIndices ]
        faceNormalCoor = self._get_face_normals(verts)
        # note: following list are the face indices that each vertex is a member.
        #       Since the number of faces a vertex is attached ranges from 1 - 6,
        #       needs a list, not an np.array.
        vfIndicesList = [ None ] *len(self.vertexCoor)
        for i in range(len(vfIndicesList)) : vfIndicesList[i] = []

        for fInx in range( len(self.fvIndices) ) :
            face = self.fvIndices[fInx]
            for vInx in face :
                vfIndicesList[vInx].append(fInx)

        # DevNote: probably a 'cleaner' method?  MROI
        normals = [ None ] *len(self.vertexCoor)
        for i in range(len(normals)) :
            vert = vfIndicesList[i]
            vertList = faceNormalCoor[ vert ]
            vertSum = np.sum(vertList,axis=0)/len(vertList)
            unitVector = np.divide( vertSum, np.linalg.norm(vertSum) )
            normals[i] = unitVector

        normals = np.array(normals)
        return normals

    def _get_face_normals(self,faceCoor,ausf=None) :
        """ Unit normals of triangular face coordinates."""

        vfT = np.transpose(faceCoor, (1,2,0) )
        vAB = np.subtract( vfT[1], vfT[0]  )
        vAC = np.subtract( vfT[2], vfT[0]  )
        vABt = np.transpose(vAB)
        vACt = np.transpose(vAC)
        if ausf is not None:
            vABt = vABt*ausf
            vACt = vACt*ausf
        cross = np.cross(vABt,vACt)
        sz = np.linalg.norm(cross,axis=1)[:,np.newaxis]
        return np.divide(cross,sz)

    def _get_face_centers(self,faceCoor) :
        """ Face centers from face vertex coordinates. """
        vfT = np.transpose(faceCoor, (1,2,0) )
        numVerts = vfT.shape[0]
        sumV = np.sum(vfT, axis=0)
        return np.transpose(sumV)/numVerts

    def _set_geometric_bounds(self) :
        """ set the bounds dictionary from surface vertex coordinates. """
        _set_geomBounds(self.vertexCoor,self._bounds)
        return    

    def _viewportCoor(self,xyzCoor,viewport=None) :
        """ UV coordinates for image texture mapping """
        # better way to do this? MROI
        xyz = np.transpose(xyzCoor)
        ab = self._map3Dto2Dsquare( xyz )
        a,b = ab
        if viewport is None : return a,b,np.full(len(a),True)
        As,Bs,Ae,Be = viewport
        Aref = 1.0-As
        Bref = 1.0-Bs
        Adelta = np.mod(Ae+Aref,1.0)
        Bdelta = np.mod(Be+Bref,1.0)
        if Adelta <=0.0 : Adelta = 1.0
        if Bdelta <=0.0 : Bdelta = 1.0
        Ap = np.mod(a+Aref,1.0)
        Bp = np.mod(b+Bref,1.0)
        Av = np.divide(Ap,Adelta)
        Bv = np.divide(Bp,Bdelta)
        inViewport = Av<=1
        inViewport = np.where(Bv>1,np.full(len(inViewport),False),inViewport)
        Av = np.clip(Av,0.0,1.0)
        Bv = np.clip(Bv,0.0,1.0)
        return Av,Bv,inViewport

    def _triangulate_faceIndices(self,fvIndx) :
        fvIndices, _ = self._triangulate_faceIndices_faceColors(fvIndx)
        return fvIndices

    def _triangulate_faceIndices_faceColors(self,fvIndx,faceColors=None) :
        """
        Subdivide all faces into triangles.  Two uses with number face edges is 3,4,5 or 6:
         1. initial triangulate collection of faces having a differnt number of edges
            called from __init__ with faceIndices passed as a list. 
         2. triangulate collection of faces with identical number of edges
            called from triangulate where faceIndicess is a numpy array
        """
        # .............................................................
        def divideFace4(face) :
            dist_a = np.linalg.norm( self.vertexCoor[ face[0]] - self.vertexCoor[ face[2]] )
            dist_b = np.linalg.norm( self.vertexCoor[ face[1]] - self.vertexCoor[ face[3]] )
            if dist_a < dist_b :
                face1 = [ face[0], face[1], face[2] ]
                face2 = [ face[2], face[3], face[0] ]
            else :
                face1 = [ face[1], face[2], face[3] ]
                face2 = [ face[3], face[0], face[1] ]               
            return face1, face2
        def divideFace5(face) :
            diag = lambda i : self.vertexCoor[ face[i]] - self.vertexCoor[ face[(i+2)%5] ]
            fc3 = lambda i : [ face[i], face[(i+1)%5], face[(i+2)%5] ]
            fc4 = lambda i : [ face[(i+2)%5], face[(i+3)%5], face[(i+4)%5], face[(i+5)%5] ]
            dist = [None]*5
            for i in range(5) : dist[i] = np.linalg.norm( diag(i) )
            minDiag = 0
            for i in range(1,5) :
                if dist[i] < dist[minDiag] : minDiag = i
            face1 = fc3(minDiag)
            otherFace = fc4(minDiag)
            face2, face3 = divideFace4(otherFace)
            return face1, face2, face3
        def divideFaces6(face) :
            # .. not the best, but easiest (could calc midpoint). MROI
            dist = [None]*3
            dist[0] = np.linalg.norm( self.vertexCoor[ face[0]] - self.vertexCoor[ face[3]] )
            dist[1] = np.linalg.norm( self.vertexCoor[ face[1]] - self.vertexCoor[ face[4]] )
            dist[2] = np.linalg.norm( self.vertexCoor[ face[2]] - self.vertexCoor[ face[5]] )
            iMin = 0
            if dist[1] < dist[0]    : iMin = 1
            if dist[2] < dist[iMin] : iMin = 2
            faceA = [ face[iMin],   face[ iMin+1 ],   face[ iMin+2 ],   face[ iMin+3 ] ]
            faceB = [ face[iMin+3], face[(iMin+4)%6], face[(iMin+5)%6], face[(iMin+6)%6] ]
            face1, face2 = divideFace4(faceA)
            face3, face4 = divideFace4(faceB)
            return face1, face2, face3, face4
        # .............................................................
        # FutDev: easy to understand but need to cleanup, MROI
        subFaceIndices = []
        fColor = []
        setColors = False
        if faceColors is not None:
            setColors = True
            if len(fvIndx) != len(faceColors) :
                # if faceColors, then this should not happen, but report if so....
                warnings.warn('Internal ERROR: tfc001')
                return None
                
        for i, face in enumerate(fvIndx) :
            numVerts = len(face)
            if numVerts == 3 :
                subFaceIndices.append(face)
                if setColors : fColor += [faceColors[i]]*1
            elif numVerts == 4 :
                f1, f2 = divideFace4(face)
                subFaceIndices.append(f1)
                subFaceIndices.append(f2)
                if setColors : fColor += [faceColors[i]]*2
            elif numVerts == 5 :
                f1, f2, f3 = divideFace5(face)
                subFaceIndices.append(f1)
                subFaceIndices.append(f2)
                subFaceIndices.append(f3)
                if setColors : fColor += [faceColors[i]]*3
            elif numVerts == 6 :
                f1, f2, f3, f4 = divideFaces6(face)
                subFaceIndices.append(f1)
                subFaceIndices.append(f2)
                subFaceIndices.append(f3)
                subFaceIndices.append(f4)
                if setColors : fColor += [faceColors[i]]*4
        
        return np.array(subFaceIndices), fColor

    def _get_edges_from_faces(self,fvIndx) :
        """
        Used in __init__ if edgeIndices=None to create edgeIndices.
        """
        # .............................................................
        def edges_from_faceList( fvIndx ) :
            # NOTE: very poor execution time. (~2 orders-of-magnitude slower)
            # face index list of lists can't be used as a Numpy array
            # since edge list of multi-type faces.
            warnings.warn('multi-type polyhedron faces: use of non-optimized edge calc.')
            keyList, edgeList = [],[]
            strkey = lambda a,b : str(max(a,b)) + '_' + str(min(a,b))
            for face in fvIndx :
                numEdges = len(face)
                for i in range(numEdges) :
                    a = face[i]
                    b = face[ (i+1)%numEdges]
                    key = strkey(a,b)
                    if key not in keyList : 
                        keyList.append(key)
                        edgeList.append([a,b])
            return np.array(edgeList)
        # .............................................................
        # following line only needed when called for mixed polygon faces
        if fvIndx.ndim == 1 : return edges_from_faceList(fvIndx)

        numEdges = fvIndx.shape[1]
        indicies = [ (i+1)%numEdges for i in range(numEdges) ]
        ofs_fvIndx = fvIndx[:,indicies]
        stepA = np.stack( (fvIndx,ofs_fvIndx), axis=-1  )
        stepB = np.reshape( stepA, (-1,2 ) )
        stepC = np.sort(stepB,axis=1)
        u, indices = np.unique(stepC, axis=0, return_index=True)
        edgeList = stepC[indices]
        return edgeList

    @property
    def _dfacecolors(self) :
        """
        Direct access to the private property _facecolors.

        Appears that the self.get_facecolor() needs a call to the
        self.do_3d_projection() which will 'sort the 2D version by view depth'.
        for render construction.  Need the original order of colors.
        BANDAID fix: update for v_1.2.0

        FutDev: define a separate private property to hold these values
                and not use the Matplotlib inherited value?   MROI
        """
        return self._facecolors

    @property
    def _dedgecolors(self) :
        """
        Direct access to the private property _edgecolor.
        """
        return self._edgecolors

    @property
    def vlim(self) :
        '''Range of values associated with color'''
        return self._bounds['vlim']

    @vlim.setter
    def vlim(self,val) :
        self._bounds['vlim'] = val
        return

    @property
    def name(self) :
        '''Descriptive identifier for the surface geometry.'''
        if self._geomName is None : return ''
        return self._geomName

    @name.setter
    def name(self, val) :
        self._geomName = val
        self.set_label( '' if val is None else val )
        return

    @property
    def cname(self) :
        '''Descriptive identifier for values indicated by color.'''
        if self.valuesName is None : return ''
        return self.valuesName

    @cname.setter
    def cname(self, val) :
        self.valuesName = val
        return

    @property
    def bounds(self) :
        """
        Dictionary of surface geometric and value ranges.
        
        Each dictionary value is a 2 float array of minimum
        and maximum values of the surface.  Keys are:

        'xlim' : x-coordinate

        'ylim' : y-coordinate

        'zlim' : z-coordinate

        'r_xy' : radial distance from the z axis

        'rorg' : radial distance from the origin

        'vlim' : value functional assignments.
            
        Values are assigned from the geometry and color mapping methods,
        including clipping.

        """

        return self._bounds

    @property
    def cBar_ScalarMappable(self) :
        """matplotlib.cm.ScalarMappable object for surface values."""
        
        vmin, vmax = self._bounds['vlim']
        norm = colors.Normalize(vmin=vmin,vmax=vmax)
        objMap = self.get_cmap()
        sm = cm.ScalarMappable(cmap=objMap, norm=norm)
        sm.set_array([])
        return sm        
    
    @property
    def edges(self) :
        """ColorLine3DCollection of the surface edges from facecolors."""
        # ......................................................
        def get_edge_color() :
            efcList = copy.copy(self.efIndicesList)
            # an edge may have 1 or 2 adjacent faces....
            for i,cList in enumerate(self.efIndicesList) :
                if len(cList)==1 : efcList[i].append(efcList[i][0])
            efcArray = np.array(efcList)
            self._postProc_surfaceColors(True)
            colorA = self.facecolors[efcArray.T[0]]
            colorB = self.facecolors[efcArray.T[1]]
            # average of adjacent face RGB colors
            return ( colorA+colorB)/2
        # ......................................................
        v = self.vertexCoor
        e = self._get_edges_from_faces(self.fvIndices)  # accounts for clipping...
        lcol = ColorLine3DCollection(v,e,'edges')
        lcol._ledgecolors = get_edge_color()  # update for v_1.2.0
        return lcol

    @property
    def initedges(self) :
        ''' ColorLine3DCollection of the initial surface edges.
            ONLY available for surfaces where the number of
            face edges differ among faces and the number of
            faces is less than 250.

        '''

        if self._initedges is None :
            raise ValueError('initedges property not available, all faces have equal number of edges (use edges property)')       
        e = self._initedges
        v = self._initverts
        lcol = ColorLine3DCollection(v,e,'edges')
        return lcol

    @property
    def vertices(self) :
        """A 3 x N array of surface vertices."""
        v = self.vertexCoor
        m = np.transpose(v)
        return m

    @property
    def facecenters(self) :
        """A 3 x N array of surface face center coordinates."""
        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        m = np.transpose(faceCenterCoor)
        return m

    @property
    def facecolors(self) :
        """A N x 4 array of surface face colors."""
        self._postProc_surfaceColors()
        return self._dfacecolors  # update for v_1.2.0

    @property
    def svd_dict(self) :
        """
        A dictionary of results from a Singular Value Decomposition of the data.
        
        Keys are:

        'disarr' : array of N floats, normalized to the surface size.
        
        'sigma' : standard deviation

        'trans' : [rotation matrix, scaling, translation]
    
        Values are assigned using the map_geom_from_svd method.
        The dictionary is None if the method has not been called by the surface object.
        
        """
        return self._svd_dict

    @property
    def area_h2b(self) :
        """ A 2 x N array of normalized N face areas and shapes."""
        f,v = self.fvIndices, self.vertexCoor
        if f.shape[1] == 3 : # triangular
            A = v[ f[:,1] ] - v[ f[:,0] ]
            B = v[ f[:,2] ] - v[ f[:,0] ]
            C = v[ f[:,2] ] - v[ f[:,1] ]
            area = 0.5*np.linalg.norm( np.cross(A,B),axis=1 )
            aveArea = np.average(area)
            normArea = area/aveArea
            edges = np.linalg.norm( np.array( [A,B,C]),axis=2)
            maxsizes = np.amax(edges,axis=0)
            f = 4*np.sqrt(3)/3 # equilateral triangle
            skew = f*area/(maxsizes*maxsizes)
        else : # quadrateral is two triangles.
            b02 = v[ f[:,2] ] - v[ f[:,0] ] 
            b13 = v[ f[:,3] ] - v[ f[:,1] ]
            base = np.array( [b02,b13] )
            diag = np.linalg.norm( base ,axis=2)
            i = np.where(diag[0]>diag[1], np.zeros(len(f),int), np.ones(len(f),int) )
            maxsizes =  diag[i,1]
            A = b02
            # top triangle ....
            B = v[ f[:,3] ] - v[ f[:,0] ]
            areaTop = 0.5*np.linalg.norm( np.cross(A,B),axis=1 )
            # bottom triangle ....
            B = v[ f[:,1] ] - v[ f[:,2] ]
            areaBtm = 0.5*np.linalg.norm( np.cross(-A,B),axis=1 )
            area = np.add( areaTop,areaBtm)
            aveArea = np.average(area)
            normArea = area/aveArea
            f = 1 # square
            skew = f*area/(maxsizes*maxsizes)
        # ....................
        return [normArea,skew]
    
    def triangulate(self, rez=0) :
        """
        PLANAR subdivision of each face into triangular faces.
        
        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the base faces
            into triangular faces, four faces per rez.
            Rez values range from 0 to 7.
            For surfaces with 4 vertices per face, each face is first
            subdivided into two faces, then recursion proceeds.
            For surfaces with 5 vertices per face, each face is first
            subdivided into three faces, then recursion proceeds.
            
        Returns
        -------
        self : surface object
        

        """
        # .............................................................
        def midVFun(vectA, vectB) :
            mid = np.add(vectA,vectB)
            mid = np.multiply(0.5,mid)
            return mid       
        # .............................................................
        
        self._postProc_surfaceColors(True)


        if self.fvIndices.shape[1] > 3 :
            # if needed, break faces into triangles before proceeding....
            subDiv = self.fvIndices.shape[1] - 2
            self.fvIndices = self._triangulate_faceIndices(self.fvIndices)

            newFaceColors = []
            for fColor in self._dfacecolors :  # update for v_1.2.0
                newFaceColors += [fColor]*subDiv
            #self.set_facecolor( np.array(newFaceColors) )
            self.set_color( np.array(newFaceColors) )  # update for v_1.3.0

        # .... update for v_1.3.0
        #indexObj ,vertexCoor, vertexValues, phongNorms = _triangulateBase(self,rez,midVFun)
        indexObj ,vertexCoor, vertexValues, phongNorms = self._triangulateBase(rez,midVFun)
        self.vertexVals = vertexValues
        self.phongNorms = phongNorms
        # .......................
        self.fvIndices = indexObj['face']
        self.vertexCoor = vertexCoor
        self.vertexVals = vertexValues
        self.evIndices = self._get_edges_from_faces(self.fvIndices)
        self.vfIndicesList = self._vertexCommonFaceIndices()
        self.efIndicesList = self._edgeCommonFaceIndices()

        v = self.vertexCoor
        verts = v[self.fvIndices]
        self.set_verts( verts )
        self._set_geometric_bounds()

        if rez != 0 :
            newFaceColors = []
            subDiv = 4**rez
            for fColor in self._dfacecolors :  # update for v_1.2.0
                newFaceColors += [fColor]*subDiv
            #self.set_facecolor( np.array(newFaceColors) )
            self.set_color( np.array(newFaceColors) )  # update for v_1.3.0

        return self        

    def normalize_scale(self) :
        """Scaling normalization array and reciprocal."""
        bounds = self.bounds
        scX = ( bounds['xlim'][1] - bounds['xlim'][0]  )
        scY = ( bounds['ylim'][1] - bounds['ylim'][0]  )
        scZ = ( bounds['zlim'][1] - bounds['zlim'][0]  )
        recip = [scX,scY,scZ]
        scale = np.reciprocal(recip)
        return scale, recip

    def vertexnormals(self, **kargs) :
        """
        Unit directions or Vector3DCollection of vertex normals at vertices.

        Parameters
        ----------
        scale: number, optional
            If not spectified, scaled proportional to the mean edge
            length of surface base faces.

        v3d: boolean, optional, default: True
            If True, Vector3DCollection is returned, else a N x 3 array
            of unit vector vertex normals is returned if False.

        color : str or float list of length 3 or 4.
            RGB or RGBA color, either a Matplotlib format string or a 
            list of color values in range [0,1].
        
        width : number, optional, default: 1
            Line width of the vector.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Returns
        -------
        Vector3DCollection object or array of unit direction vectors

        """

        kargDft = { 'scale':self._normal_scale, 'v3d':True, 'color':None, 'width':1, 'alr':_DFT_ALR  }
        KW = KWprocessor(kargDft,'vectorfield_from_op',**kargs)
        scale = KW.getVal('scale')
        isV3D = KW.getVal('v3d')
        color = KW.getVal('color')
        width = KW.getVal('width')
        alr = KW.getVal('alr')

        verts = self.vertexCoor
        norms = self._get_vertex_normals()
        if not isV3D : return norms
        name = 'vertex_normals'
        lcol = Vector3DCollection(verts,norms*scale,alr,name,colors=color, linewidths=width)
        lcol.coorType = self.coorType
        lcol.uvwOrientationType = _COORSYS["XYZ"]
        if color is None : lcol._set_vectColor( self.vertexColor )
        return lcol

    def facenormals(self, **kargs) :
        """
        Unit directions or Vector3DCollection of face normals at face centers.

        Parameters
        ----------
        scale: number, optional
            If not spectified, scaled proportional to the mean edge
            length of surface base faces.

        v3d: boolean, optional, default: True
            If True, Vector3DCollection is returned, else a N x 3 array
            of unit vector face normals is returned if False.

        color : str or float list of length 3 or 4.
            RGB or RGBA color, either a Matplotlib format string or a 
            list of color values in range [0,1].
        
        width : number, optional, default: 1
            Line width of the vector.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Returns
        -------
        Vector3DCollection object or array of unit direction vectors

        """
        #Note: default color:None so that can use as a flag to use face colors if not defined. 
        kargDft = { 'scale':self._normal_scale, 'v3d':True, 'color':None, 'width':1, 'alr':_DFT_ALR  }
        KW = KWprocessor(kargDft,'vectorfield_from_op',**kargs)
        scale = KW.getVal('scale')
        isV3D = KW.getVal('v3d')
        color = KW.getVal('color')
        width = KW.getVal('width')
        alr = KW.getVal('alr')

        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        faceNormalCoor = self._get_face_normals(verts)
        if not isV3D : return faceNormalCoor
        name = 'face_normals'
        lcol = Vector3DCollection(faceCenterCoor,faceNormalCoor*scale,alr,name,colors=color, linewidths=width)
        # --- 1.1 mod --------------------
        lcol.coorType = self.coorType
        lcol.uvwOrientationType = _COORSYS["XYZ"]
        if color is None : lcol._set_vectColor(self._dfacecolors)  # update for v_1.2.0
        return lcol

    def dispfield_from_op(self, operation, **kargs) :
        """
        Vector3DCollection of displacement of vertices to a different position.

        Parameters
        ----------
        operation : function object
            Function that takes one coordinate
            argument, a 3xN Numpy array.
            The function returns a 3xN array of vectors.
        
        returnxyz : bool { True, False }, optional, default: False
            By default, native coordinates are returned by the
            operation function.  If set True, the operation
            returns xyz Cartesian coordinates.

        useBase : bool { True, False }, optional, default: False
            When set False, vertices of the surface, prior to
            any geometric mapping or transforms, are passed to
            the operation function.  Otherwise, when True, the
            current surface vertices are passed.

        scale: number, optional, default: 1.

        color : str or float list of length 3 or 4.
            RGB or RGBA color, either a Matplotlib format string or a 
            list of color values in range [0,1].
        
        width : number, optional, default: 1
            Line width of the vector.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Returns
        -------
        Vector3DCollection object

        """

        kargDft = { 'returnxyz':False, 'useBase':False, 
                    'scale':1, 'color':'black', 'width':1, 'alr':_DFT_ALR  }
        KW = KWprocessor(kargDft,'dispfield_from_op',**kargs)
        returnxyz = KW.getVal('returnxyz')
        useBase = KW.getVal('useBase')
        scale = KW.getVal('scale')
        color = KW.getVal('color')
        width = KW.getVal('width')
        alr = KW.getVal('alr')

        if self._isCompositeSurface :
            warnings.warn('Vector operation not available for combined shapes.')
            return
        start_coor = self.vertexCoor
        if useBase : start_coor = self.baseVertexCoor
        v = self._transformVector(start_coor, operation, returnxyz)
        delta = v - self.vertexCoor
        delta = scale*delta
        # --- 1.1 mod --------------------
        vField = self._get_vectorfield_Line3DCollection(delta, color, width, alr)
        if returnxyz : vField.uvwOrientationType = _COORSYS["XYZ"]
        return vField

    def vectorfield_from_op(self, operation, **kargs) :
        """
        Vector3DCollection of vectors in u,v,w coordinates at surface vertices.

        Parameters
        ----------
        operation : function object
            Function that takes one coordinate
            argument, a 3xN Numpy array of native coordinates.
            The function returns a 3xN array of vectors.
        
        scale: number, optional, default: 1.

        color : str or float list of length 3 or 4.
            RGB or RGBA color, either a Matplotlib format string or a 
            list of color values in range [0,1].
        
        width : number, optional, default: 1
            Line width of the vector.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Returns
        -------
        Vector3DCollection object

        """

        kargDft = { 'scale':1, 'color':'black', 'width':1, 'alr':_DFT_ALR  }
        KW = KWprocessor(kargDft,'vectorfield_from_op',**kargs)
        scale = KW.getVal('scale')
        color = KW.getVal('color')
        width = KW.getVal('width')
        alr = KW.getVal('alr')

        if self._isCompositeSurface :
            warnings.warn('Vector operation not available for combined shapes.')
            return

        xyz = np.transpose(self.vertexCoor)
        abc = self.coor_convert(xyz)
        abc = np.array(abc)
        rst = np.array(operation(abc))
        rst = scale*rst
        delta = self._uvw_to_xyz(rst,abc)
        # --- 1.1 mod --------------------
        delta = np.array(delta)
        return self._get_vectorfield_Line3DCollection(delta, color, width, alr)

    def vectorfield_to_surface(self, surface, **kargs) :
        """
        Vector3DCollection of vectors from surface to surface vertex coordinates.

        Parameters
        ----------
        surface : surface object
            Surface that matches the calling surface vectors.
            Should be the same basetype and rez for surface subclasses.

        scale: number, optional, default: 1.

        color : str or float list of length 3 or 4.
            RGB or RGBA color, either a Matplotlib format string or a 
            list of color values in range [0,1].
        
        width : number, optional, default: 1
            Line width of the vector.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        Raises
        ------
        ValueError
            Mismatched surfaces based on the number of vertices.

        Returns
        -------
        Vector3DCollection object

        """

        kargDft = { 'scale':1, 'color':'black', 'width':1, 'alr':_DFT_ALR  }
        KW = KWprocessor(kargDft,'vectorfield_to_surface',**kargs)
        scale = KW.getVal('scale')
        color = KW.getVal('color')
        width = KW.getVal('width')
        alr = KW.getVal('alr')

        surLen, selfLen = len(surface.vertexCoor) , len(self.vertexCoor)
        #if surLen is not selfLen :  # 1.1 err. corrected
        if surLen != selfLen :
            raise ValueError('Surfaces have unequal number of vertices, {} != {}'.format(surLen, selfLen))            
        delta =  surface.vertexCoor - self.vertexCoor
        delta = scale*delta
        # --- 1.1 mod --------------------
        vField = self._get_vectorfield_Line3DCollection(delta, color, width, alr)
        vField.uvwOrientationType = _COORSYS["XYZ"]
        return vField

    def _map_color_from_image_vert(self, img, viewport ) :

        a,b,inViewport = self._viewportCoor(self.vertexCoor,viewport)
        
        height, width = np.subtract(img.shape[:2],1)
        M_index = ( (1-b)*height ).astype(int)
        N_index = ( a*width ).astype(int)
        #'''

        orig_colors = self.vertexColor

        if len(orig_colors) == 1 :
            onearr = np.ones( len(self.vertexCoor) )[:, np.newaxis]
            orig_colors = orig_colors*onearr       

        colorMap = []
        oneAlp = np.array([1.0])
        if viewport is None : colorMap = img[M_index,N_index]
        else :
            # FutDev: use numpy methods to optimize efficiency?
            for i in range( len(orig_colors) ) :
                colorMap.append(orig_colors[i])
                if inViewport[i] :
                    temp =  img[M_index[i],N_index[i]]
                    #.. note: img may or may not have an alpha channel..
                    temp2 = temp[0:3]
                    colorMap[i] = np.concatenate((temp2,oneAlp))   

        self.vertexColor = np.array(colorMap)
        return

    def map_color_from_image(self, fname, viewport=None) :
        """
        Assign image color values to surface face colors.

        This method is not available for composite surfaces.
        The entire image domain is applied to the coloring operation.

        Parameters
        ----------
        fname : str or file-like
            The image file to read: a filename, a URL or a file-like object opened
            in read-binary mode.  Currently restricted to PNG format.
        
        viewport : 4D array-like, optional, default: None
            Viewport defines the subdomain of the surface onto which
            the image is mapped.  The subdomain is set by normalized
            native surface coordinates in a 4D array.
            The entire surface is colored for the default of None.

        Returns
        -------
        self : surface object

        """
        
        if self._isCompositeSurface :
            warnings.warn('Image operation not available for combined shapes.')
            return self
        if self._surfaceShapeChanged :
            warnings.warn('Image operation may be anomalous after shape modification.')

        if viewport is not None :
            viewport = np.array(viewport)
            if np.any( (viewport<0) | (viewport>1) ) :
                raise ValueError('viewport list values, {}, must be between 0 and 1'.format(viewport))
            if len(viewport) != 4 :
                raise ValueError('viewport list must be length 4, found {}'.format(len(viewport)))      
        verts = self.vertexCoor[ self.fvIndices ]
        
        faceCenterCoor = self._get_face_centers(verts)
        a,b,inViewport = self._viewportCoor(faceCenterCoor,viewport)
        img = image.imread(fname)
        
        height, width = np.subtract(img.shape[:2],1)
        M_index = ( (1-b)*height ).astype(int)
        N_index = ( a*width ).astype(int)
        
        orig_colors = self._dfacecolors  # update for v_1.2.0

        if len(orig_colors) == 1 :
            onearr = np.ones( len(faceCenterCoor) )[:, np.newaxis]
            orig_colors = orig_colors*onearr       

        colorMap = []
        if viewport is None : colorMap = img[M_index,N_index]
        else :
            # FutDev: use numpy methods to optimize efficiency?
            for i in range( len(orig_colors) ) :
                colorMap.append(orig_colors[i])
                if inViewport[i] : colorMap[i] = img[M_index[i],N_index[i]]    
        self.set_color(colorMap)
        self.isSurfaceColored = True
        self.vertexValues = None
        # ========================================================================
        self._map_color_from_image_vert(img,viewport )
        # ========================================================================

        return self

    def _map_color_from_op_vert(self, operation, rgb=True) :
        xyz = np.transpose(self.vertexCoor)
        abc = self.coor_convert(xyz)
        colors = np.array(operation(abc))
        colors = np.transpose(colors)
        colors = np.clip(colors,0,1)
        if rgb : RGB = colors
        else :   RGB = cm.colors.hsv_to_rgb(colors)

        if RGB.shape[1] == 3 :  # add alpha channel.
            ones = np.array( [np.ones( RGB.shape[0] )] ).T
            RGB = np.concatenate( (RGB, ones ), axis=1  )

        self.vertexColor = RGB 
        self._bounds['vertvlim'] = _DFT_VERTVLIM

        return 

    def map_color_from_op(self, operation, rgb=True, cname=None) :
        """
        Assignment of face color from a function.

        Face colors are assigned from a function
        of face coordinates.

        Parameters
        ----------
        operation : function object
            Function that takes one argument,
            a 3xN Numpy array of native coordinates.
            The function returns a 3xN color value.

        rgb : bool {True, False}, optional, default: True
            By default, RGB color values are returned by the
            operation function.  If set False, the operation
            returns HSV color values.

        cname : string, optional, default: None or op function name.
            Descriptive identifier for values indicated by color.

        Returns
        -------
        self : surface object

        """

        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        xyz = np.transpose(faceCenterCoor)
        abc = self.coor_convert(xyz)
        colors = np.array(operation(abc))
        colors = np.transpose(colors)
        colors = np.clip(colors,0,1)
        if rgb : RGB = colors
        else :   RGB = cm.colors.hsv_to_rgb(colors)
        self.set_color(RGB)
        self.isSurfaceColored = True
        cname = _getFunctionName(operation,cname) 
        if cname is not None : self.valuesName = cname 
        self._bounds['vlim'] = _DFT_VLIM

        self._map_color_from_op_vert(operation, rgb)

        return self

    def _map_cmap_from_datagrid_vert(self, datagrid, cmap=None, viewport=None) :

        a,b,inViewport = self._viewportCoor(self.vertexCoor,viewport)      
        data = datagrid
        dmax = np.amax(datagrid)
        dmin = np.amin(datagrid)
        xd = np.linspace(0, 1, data.shape[0] )
        yd = np.linspace(0, 1, data.shape[1] )
        g =  interpolate.interp2d(yd, xd, data, kind='cubic')
        d = []
        # FutDev: use numpy methods to optimize efficiency.
        for ab in np.transpose([a,b]) :
            d.append( g(ab[0],ab[1])[0] )   
        d = np.where(inViewport,d,np.full(len(d),dmin))
        
        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cmap = self.get_cmap()
        
        orig_colors = self.vertexColor

        if len(orig_colors) == 1 :
            onearr = np.ones( len(self.vertexCoor) )[:, np.newaxis]
            orig_colors = orig_colors*onearr       
        
        norm = colors.Normalize(dmin,dmax)
        temp = cmap(norm(d))
        if viewport is None :
            colorMap = temp
        else :
            colorMap = np.where(inViewport[:, np.newaxis], temp, orig_colors )

        # FutDev: set_color used instead of set_facecolor to eliminate 'gaps' between faces,  <<<<<
        # however, set edgecolor required after illumination if edges are to be displayed.    <<<<<
        # HOWEVER !!!!  if colorMap has a alpha<1, these don't register with the edges.       <<<<<

        #self.vertexValues = d   <<<< NO.. interpolated values only in viewport,
        self.vertexColor = np.array(colorMap)

        return

    def map_cmap_from_datagrid(self, datagrid, cmap=None, viewport=None, cname=None) :
        """
        Face color assignment using a 2D datagrid surface.
        
        Datagrid values are normalized in the range 0 to 1.

        This method is not available for composite surfaces.
        The entire datagrid domain is applied to the geometric operation.

        Parameters
        ----------
        datagrid : 2D float array

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the datagrid values to colors.

        viewport : 4D array-like, optional, default: None
            Viewport defines the subdomain of the surface onto which
            the datagrid is mapped.  The subdomain is set by normalized
            native surface coordinates in a 4D array.
            The entire surface is colored for the default of None.

        cname : string, optional, default: None.
            Descriptive identifier for values indicated by color.

        Returns
        -------
        self : surface object

        """

        if self._isCompositeSurface :
            warnings.warn('Datagrid operation not available for combined shapes.')
            return self
        if self._surfaceShapeChanged :
            warnings.warn('Datagrid operation may be anomalous after shape modification.')

        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        a,b,inViewport = self._viewportCoor(faceCenterCoor,viewport)      
        data = datagrid
        dmax = np.amax(datagrid)
        dmin = np.amin(datagrid)
        xd = np.linspace(0, 1, data.shape[0] )
        yd = np.linspace(0, 1, data.shape[1] )
        g =  interpolate.interp2d(yd, xd, data, kind='cubic')
        d = []
        # FutDev: use numpy methods to optimize efficiency.
        for ab in np.transpose([a,b]) :
            d.append( g(ab[0],ab[1])[0] )   
        d = np.where(inViewport,d,np.full(len(d),dmin))
        
        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cmap = self.get_cmap()
        
        orig_colors = self._dfacecolors  # update for v_1.2.0

        if len(orig_colors) == 1 :
            onearr = np.ones( len(faceCenterCoor) )[:, np.newaxis]
            orig_colors = orig_colors*onearr       
        
        norm = colors.Normalize(dmin,dmax)
        temp = cmap(norm(d))
        if viewport is None :
            colorMap = temp
        else :
            colorMap = np.where(inViewport[:, np.newaxis], temp, orig_colors )

        # FutDev: set_color used instead of set_facecolor to eliminate 'gaps' between faces,
        # however, set edgecolor required after illumination if edges are to be displayed.
        # HOWEVER !!!!  if colorMap has a alpha<1, these don't register with the edges.
        self.set_color(colorMap)
        self.isSurfaceColored = True
        self._bounds['vlim'] = [ dmin, dmax ]

        self.valuesName = cname
        # ========================================================================
        self._map_cmap_from_datagrid_vert(datagrid, cmap, viewport )
        # ========================================================================

        return self

    def _getColorMap_fromNormals(self,fvo,cmap, direction, isAbs) :
        """ used by map_cmap_from_normals() & _map_cmap_from_normals_vert() """
        unitVector = lambda v : np.divide( v, np.linalg.norm(v) )
        incidentLight = unitVector( direction )
        d = np.dot(fvo,incidentLight)
        if isAbs : v = 1 - np.abs(d)
        else :
            v = ( d + 1 )/2
        v = np.abs(v)  # to insure positive values near zero
        colorMap = cmap(v)
        return colorMap

    def _map_cmap_from_normals_vert(self, cmap, direction, isAbs) :
        vobj = self._get_vertex_normals()
        cMap = self._getColorMap_fromNormals(vobj, cmap, direction, isAbs)
        self.vertexColor = np.array(cMap)
        return

    def map_cmap_from_normals(self, cmap=None, direction=None, cname=None, isAbs=False) :
        """
        Face color assignment using normals relative to a direction.

        The dot product of face normals with the direction 
        is used to assign face colors from a colormap.

        Parameters
        ----------
        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the dot product values to colors.
        
        direction : list of size 3, optional, default: [1,0,1]
            A 3D vector in xyz Cartesian coordinates designating
            the direction of the illumination source.
            If assigned the Axes3D, will use the view direction.

        cname : string, optional, default: None.
            Descriptive identifier for values indicated by color.

        isAbs : boolean, optional, default: False
            If set True, the absolute value of the dot product
            between face normals and illumination direction is
            used.

        Returns
        -------
        self : surface object

        """
        if direction is not None:
            temp = direction
            try :  # check if direction holds the elev, azim (ie. 3D axis)
                elev,azim = temp.elev, temp.azim
                direction = elev_azim_2vector(elev,azim)
            except :
                direction = temp
                if len(direction) != 3 :
                    raise ValueError('Cmap from normals direction must be a array of length 3, or a Matplotlib Axes3D object.')
        else : direction = _ILLUM

        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cm_colorMap = self.get_cmap()

        verts = self.vertexCoor[ self.fvIndices ]
        fvobj = self._get_face_normals(verts)
        cMap = self._getColorMap_fromNormals(fvobj, cm_colorMap, direction, isAbs)
        self.set_cmap(cm_colorMap)
        self.set_color(cMap)
        self.isSurfaceColored = True
        self.onlyColoredFaces = True

        self.valuesName = cname
        # ========================================================================
        self._map_cmap_from_normals_vert( cm_colorMap, direction, isAbs )
        # ========================================================================

        return self

    def map_cmap_from_vertvals(self,cmap=None, cname=None) :
        """
        Face color assignment using face vertex values.

        Parameters
        ----------
        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the datagrid values to colors.
        
        cname : string, optional, default: None.
            Descriptive identifier for values indicated by color.

        Returns
        -------
        self : surface object

        """
        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cmap = self.get_cmap()

        if self.vertexVals is None :
            warnings.warn('ERROR: [map_cmap_from_vertvals] '+
                    'vertice values are unassigned, no action taken.')
            return self

        fvIndices = self.fvIndices
        vertVals = self.vertexVals
        fvIndices = np.array(fvIndices)
        fvertVals = vertVals[fvIndices]
        v = np.average(fvertVals,axis=1)
        vmin, vmax = v.min(), v.max()
        self.faceValues = v
        self._bounds['vlim'] = [ vmin, vmax ]
        norm = colors.Normalize(vmin=vmin,vmax=vmax)
        colorMap  = cmap(norm(v))
        self.set_color(colorMap )
        self.isSurfaceColored = True
        self.valuesName = cname
        return self

    def set_vertvals(self,values) :
        """
        Assign scalar values to each vertex.

        Parameters
        ----------
        values : list or array of length N, where N is the
            number of vertices.

        Returns
        -------
        self : surface object

        """
        values = np.array(values)
        noVerts = self.vertexCoor.shape[0]
        noVals  = values.shape[0]
        if noVals != noVerts :
            warnings.warn("Error: [set_vertvals] no. vertices {} != no. values {}, no action taken".format(noVerts,noVals))
            return self
        self.vertexVals = values
        return self

    def _map_cmap_from_op_vert(self, operation, cmap) :
        
        xyz = np.transpose(self.vertexCoor)
        # .....
        abc = self.coor_convert(xyz)
        v = np.array(operation(abc))
        self.vertexValues = v
        norm = colors.Normalize(vmin=v.min(),vmax=v.max())
        colorMap = cmap(norm(v))
        self.vertexColor = colorMap
        colorMap = np.array(colorMap)

        self._bounds['vertvlim'] = [ v.min(), v.max() ]

        return 

    def map_cmap_from_op(self, operation=None, cmap=None, cname=None) :
        """
        Functional assignment of a color from a color map.

        Face coordinates are used to calculate a scalar 
        which is then used to assign face colors from a
        colormap.

        Parameters
        ----------
        operation : function object, default : None
            Function that takes one argument,
            a 3xN Numpy array of native coordinates.
            The function returns a Numpy array of scalar values.
            If function is None, function will map in the z-direction
            or r-direction, dependent on the native coordinates.

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the function return values to colors.
        
        cname : string, optional, default: None or op function name.
            Descriptive identifier for values indicated by color.

        Returns
        -------
        self : surface object

        """

        dftDirName = None
        if operation is None :
            opDir = 2
            dftDirName = 'Z-direction'
            if self.coorType == _COORSYS["SPHERICAL"] :   opDir = 0
            if self.coorType == _COORSYS["CYLINDRICAL"] : opDir = 0
            if opDir == 0 : dftDirName = 'R-direction'
            operation = lambda c : c[opDir]

        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cmap = self.get_cmap()

        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)
        xyz = np.transpose(faceCenterCoor)
        # .....
        abc = self.coor_convert(xyz)
        v = np.array(operation(abc))
        self.faceValues = v
        vmin, vmax = v.min(), v.max()
        norm = colors.Normalize(vmin=v.min(),vmax=v.max())
        
        colorMap = cmap(norm(v))
        # FutDev: set_color used instead of set_facecolor to eliminate 'gaps' between faces,
        # however, set edgecolor required after illumination if edges are to be displayed.
        # HOWEVER !!!!  if colorMap has a alpha<1, these don't register with the edges.
        self.set_color(colorMap)
        self.isSurfaceColored = True
        self._bounds['vlim'] = [ vmin, vmax ] 

        #cname = _getFunctionName(operation,cname)
        trialName = None
        if cname is None :
            if dftDirName is not None : trialName = dftDirName
            else :  trialName =  _getFunctionName(operation,None)
        else: trialName = cname
        self.valuesName = trialName
        # ========================================================================
        self._map_cmap_from_op_vert(operation, cmap)
        # ========================================================================
        
        return self

    def map_geom_from_datagrid(self,datagrid, scale=1.0, viewport=None, name=None ) :
        """
        Append surface vertices proportional to 2D datagrid surface.
        
        Datagrid values are normalized in the range 0 to 1.

        This method is not available for composite surfaces.
        The entire datagrid domain is applied to the geometric operation.

        Parameters
        ----------
        datagrid : 2D array

        scale : number, optional, default: 1

        viewport : 4D array-like, optional, default: None
            Viewport defines the subdomain of the surface onto which
            the datagrid is mapped.  The subdomain is set by normalized
            native surface coordinates in a 4D array.
            The entire surface is geometrically mapped for the default of None.
        
        name : string, optional, default: None
            Descriptive identifier for the geometry.

        Returns
        -------
        self : surface object

        """

        if self._isCompositeSurface :
            warnings.warn('Datagrid operation not available for combined shapes.')
            return self
        if self._surfaceShapeChanged :
            #warnings.warn('Datagrid operation may be anomalous after shape modification.')
            warnings.warn('Datagrid operation not available after shape modification.')
            return self

        a,b,inViewport = self._viewportCoor(self.vertexCoor,viewport)

        data = datagrid
        dmax = np.amax(datagrid)
        dmin = np.amin(datagrid)
        delta = dmax-dmin
        data = (datagrid - dmin)/delta
        xd = np.linspace(0, 1, data.shape[0] )
        yd = np.linspace(0, 1, data.shape[1] )
        g =  interpolate.interp2d(yd, xd, data, kind='cubic')

        # FutDev: optimize efficiency (?).
        d = []
        for ab in np.transpose([a,b]) :
            d.append( g(ab[0],ab[1])[0] )
        d = np.where(inViewport,d,np.zeros(len(d)))
        d = d[:, np.newaxis]

        # FutDev: 'fVal' only applicable for initial subclassed
        # surfaces to get surface normals at a vertex.
        # Future work: to make this generic for any surface.
        fVal = np.array([0.,0.,1.])
        if self.disp_Vector is not None :
            fVal = self.vertexCoor*self.disp_Vector

        v =  self.vertexCoor + scale*d*fVal
        verts = v[self.fvIndices]
        self.vertexCoor = v
        self._set_geometric_bounds()
        self.set_verts( verts )
        self._surfaceShapeChanged = True
        if name is not None: self._geomName = name
        self.phongNorms = None  # update for v_1.3.0
        return self

    def map_geom_from_image(self, fname, scale=1.0, viewport=None, cref='v', hzero=0, name=None ) :
        """
        Append surface vertices proportional to image color component values.

        For planar and polar coordinates, the z-coordinate is appended.
        For cylindrical and spherical coordinates, the r-coordinate is appended.

        This method is not available for composite surfaces.
        The entire image domain is applied to the geometric operation.

        Parameters
        ----------
        fname : str or file-like
            The image file to read: a filename, a URL or a file-like object opened
            in read-binary mode.  Currently restricted to PNG format.

        scale : number

        cref : {'r','g','b','h','s','v'}, optional, default: 'v'
            Sets which color tuple value to use for geometry mapping.
            Note: only the first character of the string is evaluated.

        hzero : number or None, optional, default: 0
            Argument is used if cref is set to 'h'.
            The hzero magnitude indicates the Hue value for a zero coordinate
            diplacement. The hzero sign (positive or negative),
            indicates the direction of increasing value in the
            appended coordinate with hue value.  Range of hzero is [-1,1].

        viewport : 4D array-like, optional, default: None
            Viewport defines the subdomain of the surface onto which
            the image is mapped.  The subdomain is set by normalized
            native surface coordinates in a 4D array.
            The entire surface is geometrically mapped for the default of None.

        name : string, optional, default: None
            Descriptive identifier for the geometry.

        Returns
        -------
        self : surface object

        """
        # ..............................................................
        def getIndex(cref) :
            cref = cref[0].lower()
            valIndex = 'rgbhsv'.find(cref)
            if valIndex < 0 :
                raise ValueError('Invalid cref argument: ' + str(cref))
            isHSV = False
            if valIndex > 2 :
                valIndex -= 3
                isHSV = True
            return valIndex, isHSV
        # ..............................................................
        if self._isCompositeSurface :
            warnings.warn('Image operation not available for combined shapes.')
            return self
        if self._surfaceShapeChanged :
            warnings.warn('Image operation may be anomalous after shape modification.')
        
        if (hzero<-1) or (hzero>1) :
            raise ValueError('hzero must be between -1 and 1, found {}'.format(hzero))
        ho = (1.0 -hzero)
        hsgn = np.sign(hzero)
        h_ref = lambda h : np.mod( ho + hsgn*h , 1.0)
        # ...................
        a,b,inViewport = self._viewportCoor(self.vertexCoor,viewport)
        img = image.imread(fname)
        height, width = np.subtract(img.shape[:2],1)
        x_index = ( (1-b)*height ).astype(int)
        y_index = ( a*width ).astype(int)

        valIndex, isHSV = getIndex(cref)
        rgbcolor = img[x_index,y_index]
        if isHSV :
            hsvcolor = colors.rgb_to_hsv(rgbcolor)
            hsvcolor = np.transpose(hsvcolor)
            if valIndex==0 : d= h_ref(hsvcolor[0])
            else :           d = hsvcolor[valIndex]
        else :
            rgbcolor = np.transpose(rgbcolor)
            d = rgbcolor[valIndex]

        d = np.where(inViewport,d,np.zeros(len(d)))
        d = d[:, np.newaxis]
        
        # FutDev: 'fVal' only applicable for initial subclassed surfaces
        # surfaces to get surface normals at a vertex.
        # Future work: to make this generic for any surface.
        fVal = np.array([0.,0.,1.])
        if self.disp_Vector is not None :
            fVal = self.vertexCoor*self.disp_Vector

        v =  self.vertexCoor + scale*d*fVal
        verts = v[self.fvIndices]
        self.vertexCoor = v
        self._set_geometric_bounds()
        self.set_verts( verts )
        self._surfaceShapeChanged = True
        if name is not None: self._geomName = name
        self.phongNorms = None  # update for v_1.3.0
        return self
    
    def map_geom_from_op(self, operation, returnxyz=False, name=None) :
        """
        Functional transformation of surface vertex coordinates.

        This method is not available for composite surfaces.

        Parameters
        ----------
        operation : function object
            Coordinate mapping function that takes one
            argument, a 3xN Numpy array of native coordinates.
            The function returns a 3xN array of coordinates.

        returnxyz : bool { True, False }, optional, default: False
            By default, native coordinates are returned by the
            operation function.  If set True, the operation
            returns xyz Cartesian coordinates.

        name : string, optional, default: None or op function name.
            Descriptive identifier for the geometry.

        Returns
        -------
        self : surface object

        """

        return self._map_geom_from_op( operation, returnxyz, name, check=True)

    def _map_geom_from_op(self, operation, returnxyz, name, check) :
        
        # used by domain with check=False to allow composites scaling.
        # otherwise, not allowed by map_geom_from_op method.

        if self._isCompositeSurface and check :
            warnings.warn('Map operation MUST use xyz coordinates for combined shapes.')
            #warnings.warn('Map operation not available for combined shapes.')
            #return self
        v = self._transformVector(self.vertexCoor, operation, returnxyz)
        verts = v[self.fvIndices]
        self.vertexCoor = v
        self._set_geometric_bounds()
        self.set_verts( verts )
        self._surfaceShapeChanged = True
        if name is None : name = self._geomName
        name = _getFunctionName(operation,name)
        if name is not None: self._geomName = name
        self.phongNorms = None  # update for v_1.3.0
        return self
    
    def clip(self, operation, usexyz=False) :
        """
        Remove faces from the surface based on position.

        Parameters
        ----------
        operation : function object
            Function that takes one argument, a 3xN Numpy array of coordinates.
            The function returns a bool { True, False } indicating if
            the face, at the face centered coordinate, is to be 
            retained (True), otherwise, the face is removed from the
            surface (False).

        usexyz : bool { True, False }, default: False
            If True, face centers are passed to the operation function
            in xyz coordinates.  If False, face centers are passed
            in native coordinates.

        Returns
        -------
        self : surface object

        """

        # Clip operation MUST use xyz coordinates for combined shapes.
        # since native coordinate may be undefined.
        if self._isCompositeSurface : usexyz = True

        verts = self.vertexCoor[ self.fvIndices ]
        faceCenterCoor = self._get_face_centers(verts)

        self._postProc_surfaceColors()
        orig_colors = self._dfacecolors  # update for v_1.2.0

        xyz = np.transpose(faceCenterCoor)
        if usexyz is False :  xyz = self.coor_convert(xyz)
        abc = xyz
        shouldKeep = np.array(operation(abc))
        if np.any(shouldKeep) is None :
            warnings.warn('WARNING: Clipping resulted in no faces found, surface return unclipped.')
            return self

        fvi = self.fvIndices[shouldKeep]
        fcol = orig_colors[shouldKeep]

        # this permits multiple clips without changing original, so restore is available.
        if self.restoredClip is None :
            self.restoredClip = { 'faces' : self.fvIndices, 'colors' : orig_colors }

        fvi = np.array(fvi)
        self.fvIndices = fvi
        self.set_color(fcol)        
        v = self.vertexCoor
        vfvi = v[fvi]
        self.set_verts( vfvi )

        clipVerts  = np.reshape(vfvi,(-1,3))
        _set_geomBounds( clipVerts,self.bounds )
        self._set_index_relations()
        self.phongNorms = None  # update for v_1.3.0
        return self

    def clip_alpha(self,alphaCut,useval=False) :
        """
        Remove faces from the surface based on alpha transparency.

        Parameters
        ----------
        alphaCut : scalar
            If the face color alpha is greater than alphaCut, the face is 
            retained, otherwise, the face is removed from the surface.

        useval : bool { True, False }, default: False
            If True, clipping based on color HSV value.  Otherwise,
            clipping is based on the color alpha value.

        Returns
        -------
        self : surface object

        """
        self._postProc_surfaceColors()
        orig_colors = self._dfacecolors  # update for v_1.2.0

        if useval :
            orig_hsv = colors.rgb_to_hsv(orig_colors[:,:3])
            alphas = np.ravel( orig_hsv[:,2:3] )
        else :
            alphas = np.ravel( orig_colors[:,3:4] )
        cuts = np.full(len(alphas),alphaCut)
        shouldKeep = alphas>cuts

        if np.any(shouldKeep) is None :
            warnings.warn('WARNING: Alpha clipping resulted in no faces found, surface return unclipped.')
            return self
        fvi = self.fvIndices[shouldKeep]
        fcol = orig_colors[shouldKeep]

        # this permits multiple clips without changing original, so restore is available.
        if self.restoredClip is None :
            self.restoredClip = { 'faces' : self.fvIndices, 'colors' : orig_colors }

        fvi = np.array(fvi)
        self.fvIndices = fvi
        self.set_color(fcol)        
        v = self.vertexCoor
        vfvi = v[fvi]
        self.set_verts( vfvi )

        clipVerts  = np.reshape(vfvi,(-1,3))
        _set_geomBounds( clipVerts,self.bounds )
        self._set_index_relations()
        self.phongNorms = None  # update for v_1.3.0
        return self

    def clip_plane(self,dist,**kargs) :
        """
        Remove faces from the surface based on a clip surface.

        Parameters
        ----------
        dist : float
            Distance from the origin to the clip surface,
            along the direction vector.
        
        direction : array-like, optional, default: [0,0,1]
            A xyz vector normal to the intersection plane
            for a planar clip surface or axial direction of
            a cylinder for a cylindrical clip surface.

        coor : integer or string indicating the type of clip surface:
            0, p, P, xyz,planar                 - planar (default)
            1, c, C, cylinder,pplar,cylindrical - cylinder
            2, s, S, sphere,spherical           - sphere
            3, x, X                             - y-z plane
            4, y, Y                             - x-z plane
            5, z, Z                             - x-y plane

        Returns
        -------
        self : surface object

        """
        # ...........................................................
        def clip_op(xyz,dist,coor,direction) :
            abc = np.transpose(xyz)
            dotprod = _coor_dotprod(direction,coor,abc)
            if coor > 0 :
                if dist < 0 :
                    shouldclip = np.greater(dotprod,np.full(len(dotprod),-dist))
                else :
                    shouldclip = np.less(dotprod,np.full(len(dotprod),dist))
            else :    
                shouldclip = np.less(dotprod,np.full(len(dotprod),dist))
            return shouldclip
        # ...........................................................
        dirDft = { 'direction': [0,0,1.0] }
        kargDft = { **_COOR_KWARGS, **dirDft }
        KW = KWprocessor( kargDft, 'clip_plane', **kargs)
        coor = KW.getVal('coor')
        direction = KW.getVal('direction')

        if coor==3 : coor, direction = 0, [ 1,0,0 ]
        if coor==4 : coor, direction = 0, [ 0,1,0 ]
        if coor==5 : coor, direction = 0, [ 0,0,1 ]

        return self.clip( lambda xyz : clip_op(xyz,dist,coor,direction),True)

    def clip_normals(self,direction=None) :
        """ 
        Remove faces from the surface based on a direction relative to face normals.

        Parameters
        ----------
        direction : array-like or 3Daxis, optional, default: (30,-60)
            The vector direction for viewing (elev,azim) or 3Daxis. 

        Returns
        -------
        self : surface object

        """

        if direction is not None:
            temp = direction
            try :  # check if direction holds the elev, azim (ie. 3D axis)
                elev,azim = temp.elev, temp.azim
                direction = np.array( elev_azim_2vector(elev,azim) )
            except :
                direction = temp
                if len(direction) != 3 :
                    raise ValueError('Clip from normals direction must be a array of length 3, or a Matplotlib Axes3D object.')
        else : direction = np.array( elev_azim_2vector( *_DFT_VIEW ) )

        verts = self.vertexCoor[ self.fvIndices ]
        faceNormalCoor = self._get_face_normals(verts)

        self._postProc_surfaceColors()
        orig_colors = self._dfacecolors  # update for v_1.2.0

        shouldKeep = np.dot(faceNormalCoor,direction) > 0.0

        if np.any(shouldKeep) is None :
            warnings.warn('WARNING: Normal clipping resulted in no faces found, surface return unclipped.')
            return self
        fvi = self.fvIndices[shouldKeep]
        fcol = orig_colors[shouldKeep]

        # this permits multiple clips without changing original, so restore is available.
        if self.restoredClip is None :
            self.restoredClip = { 'faces' : self.fvIndices, 'colors' : orig_colors }

        fvi = np.array(fvi)
        self.fvIndices = fvi
        self.set_color(fcol)        
        v = self.vertexCoor
        vfvi = v[fvi]
        self.set_verts( vfvi )

        clipVerts  = np.reshape(vfvi,(-1,3))
        _set_geomBounds( clipVerts,self.bounds )
        self._set_index_relations()
        self.phongNorms = None  # update for v_1.3.0
        return self

    def _restore_from_clip(self) :
        """
        Reset face and colors prior to any clip operation.
        (in development)
        """
        
        # FutDev: this method needs thorough testing.
        # probably doesn't work as is.
        if self.restoredClip is None :
            warnings.warn('There is no clipped surface to restore.')
            return self
        fvi = self.restoredClip['faces']
        fcol = self.restoredClip['colors']
        self.restoredClip = None
        self.fvIndices = fvi
        self.set_color(fcol)        
        v = self.vertexCoor
        self.set_verts(v[ np.array(fvi) ])
        self._set_geometric_bounds() # reset to original bounds...
        return self

    def set_surface_alpha(self, alpha, constant=False, adjustlw=True) :
        """
        Adjust the face color alpha values and linewidth of the surface.

        Parameters
        ----------
        alpha : scalar
              Alpha is in the range 0 to 1.

        constant : bool { True, False }, optional, False
            If False, face color values are multiplied by
            alpha.  If True, all face colors alpha channels
            are assigned to alpha.

        adjustlw : bool { True, False }, optional, True
            If True, the surface linewidth will be set
            to make the edges visually transparent with
            the faces for Matplotlib rendering.
            
        Returns
        -------
        self : surface object

        """

        if (alpha<0.0) or (alpha>1.0) :
            raise ValueError('surface alpha must be between 0 and 1, found {}'.format(alpha))

        self._postProc_surfaceColors()
        orig_colors = self._dfacecolors  # update for v_1.2.0
        if constant :
            np.put_along_axis(orig_colors, np.array([[3]]), alpha, axis=1)
            colorMap = orig_colors
        else:
            colorMap = np.multiply(orig_colors,np.array([1.0,1.0,1.0,alpha]))
        self.set_color(colorMap)      
        if adjustlw : 
            # heuristic approach to hiding face edges.  BANDAID
            lw =  0.0063*np.exp(4.2*alpha)
            self.set_linewidth(lw)
        return self   
   
    def _get_face_normals_from_phongNorms(self) :
        # called from shade and hilite  # update for v_1.3.0
        # note: lint will flag a problem : 'self.phongNorms is unsubscriptable'.
        #       self.phongNorms is set to None as a flag, but on triangulation,
        #       it is assigned an array (so no problem ;-).
        vNorms = self.phongNorms[ np.array(self.fvIndices) ]
        sumNorms = np.sum(vNorms,axis=1)
        lnorm = np.linalg.norm(sumNorms,axis=1)
        lnorm3 = np.empty([3,lnorm.shape[0]])
        lnorm3[:] = lnorm 
        norm = np.divide( sumNorms, lnorm3.T) 
        return norm

    def shade(self, depth=0, direction=None, contrast=None, isAbs=False, ax=None, rview=False, flat=True) :
        """
        Reduce surface HSV color Value based on face normals.
        
        The dot product of face normals with the illumination direction 
        is used to adjust face HSV color value.

        Parameters
        ----------
        depth : scalar, optional, default: 0
            Minimum color value of shaded surface faces.
            Depth value ranges from 0 to 1.

        direction : array-like, optional, default: (1,0,1)
            A xyz vector pointing to the illumination source. 

        contrast : scalar, optional, default: 1
            Shading contrast adjustment from low to high with a value
            of 1 for linear variations with the normal direction.
            Contrast value ranges from 0.1 to 3.     
        
        isAbs : bool, optional, default: False
            If True, the absolute value of the dot product is
            is used to determine color value.

        ax : Matplotlib 3D axes.

        rview : boolean, default: False
            If True, direction is relative to the view,
            otherwise, relative to the axes.

        flat : boolean, default: True
            If False, smooth flat triangulated faces.

        Returns
        -------
        self : surface object

        """

        if direction is None : direction = _ILLUM
        viewdir,ausf = None,None
        if ax is not None:
            ausf = _axisUSF(ax)
            viewdir = elev_azim_2vector(ax.elev,ax.azim)
            if rview : direction = rtv(direction,ax.elev,ax.azim)

        if np.any( (depth<0) | (depth>1) ) :
            raise ValueError('depth values, {}, must be between 0 and 1'.format(depth))
        if contrast is not None:
            if (contrast<0.1 or contrast>3) :
                raise ValueError('contrast must be between 0.1 and 3. , found {}'.format(contrast))
    
        # for triangulation smoothing  # update for v_1.3.0
        if not (self.phongNorms is None or flat) :
            norms = self._get_face_normals_from_phongNorms()
        else :
            verts = self.vertexCoor[ np.array(self.fvIndices) ]
            norms = self._get_face_normals(verts,ausf)

        unitDirection = np.divide( direction, np.linalg.norm(direction) )
        dprod = np.dot( norms, unitDirection )

        if viewdir is not None :
            vdot = np.dot(norms,viewdir)
            dprod = np.where(np.less(vdot,[0.0]),-dprod,dprod)

        # 0<d<1, heuristic function for 
        # normalized domain of effective dot product, d=f(dprod)
        if isAbs : d = 1 - np.abs(dprod)
        else :
            d = ( dprod + 1 )/2
        d = np.abs(d)  # to insure positive values near zero

        if contrast is not None :
            d =  0.5*( 1 + np.sign(dprod)*np.power( np.abs(dprod) , 1/contrast) )      
        
        # extract HSV, leaving alpha unchanged.
        self._postProc_surfaceColors()
        orig_colors = self._dfacecolors  # update for v_1.2.0
        alphas = orig_colors[:,3:4]
        fc_less_alpha = orig_colors[:,:3]
        hsv_vals = cm.colors.rgb_to_hsv(fc_less_alpha)

        # adjust color value with multiplier.
        cvm = ((1-depth)*d + depth*np.ones(len(d)))[:, np.newaxis]
        hsv_vals[:,2] = (cvm*hsv_vals[:,2:3])[:,0]
        rgb_vals =  cm.colors.hsv_to_rgb(hsv_vals)
        shade_colors = np.concatenate((rgb_vals, alphas), axis=1)
        self.set_color(shade_colors)  # NOTE: this will also set appropriate edge colors.
        return self

    def hilite(self, height=1, direction=None, focus=None, ax=None, rview=False, flat=True) :
        """
        Increase surface HSV color Value and reduce Saturation, based on face normals.

        The dot product of face normals with the illumination direction 
        is used to adjust face HSV color value and saturation.
        Faces having normals with a negative component in the direction of
        illumination are not affected. The reduction of color saturation is
        proportional to the increase in color value.

        Parameters
        ----------
        height : scalar, optional, default: 1
            Maximum color value of highlighted surface faces.
            Height value ranges from 0 to 1.

        direction : array-like, optional, default: (1,0,1)
            A xyz vector pointing to the illumination source. 

        focus : scalar, optional, default: 1
            Highlighting focus adjustment from low to high with a value
            of 1 for quadratic variations with the normal direction.
            Focus value ranges from 0.1 to 3.     
        
        ax : Matplotlib 3D axes.

        rview : boolean, default: False
            If True, direction is relative to the view,
             otherwise, relative to the axes.

        flat : boolean, default: True
            If False, smooth flat triangulated faces.
            
        Returns
        -------
        self : surface object

        """

        if direction is None : direction = _ILLUM
        viewdir,ausf = None,None
        if ax is not None:
            ausf = _axisUSF(ax)
            viewdir = elev_azim_2vector(ax.elev,ax.azim)
            if rview : direction = rtv(direction,ax.elev,ax.azim)

        if np.any( (height<0) | (height>1) ) :
            raise ValueError('height values, {}, must be between 0 and 1'.format(height))
        if focus is not None:
            if (focus<0.1 or focus>3) :
                raise ValueError('focus must be between 0.1 and 3. , found {}'.format(focus))
    
        # for triangulation smoothing  # update for v_1.3.0
        if not (self.phongNorms is None or flat) :
            norms = self._get_face_normals_from_phongNorms()
        else :
            verts = self.vertexCoor[ np.array(self.fvIndices) ]
            norms = self._get_face_normals(verts,ausf)

        unitDirection = np.divide( direction, np.linalg.norm(direction) )
        dprod = np.dot( norms, unitDirection )

        if viewdir is not None :
            vdot = np.dot(norms,viewdir)
            dprod = np.where(np.less(vdot,[0.0]),-dprod,dprod)

        # 0<d<1, heuristic function for 
        # normalized domain of effective dot product, d=f(dprod)
        d = np.where(dprod<0 , 0, dprod)       
        if focus is None : focus = 1
        d0 = np.power(focus,1/2)/2.0
        y = (d-d0)/(1-d0)
        d = np.where(y<0 , 0, y)
        d = np.power(d,1.0+2*focus)

        # extract HSV, leaving alpha unchanged.
        self._postProc_surfaceColors()
        orig_colors = self._dfacecolors  # update for v_1.2.0
        alphas = orig_colors[:,3:4]
        fc_less_alpha = orig_colors[:,:3]
        hsv_vals = cm.colors.rgb_to_hsv(fc_less_alpha)

        # adjust color saturation and value, linear with normalized domain.
        hd = height*d
        omhd = (np.ones( len(norms) ) - hd)
        omhd = omhd[:,np.newaxis]
        hsv_vals[:,1] = (omhd*hsv_vals[:,1:2])[:,0]
        hsv_vals[:,2] = (omhd*hsv_vals[:,2:3])[:,0]
        hsv_vals[:,2] += hd
        rgb_vals =  cm.colors.hsv_to_rgb(hsv_vals)
        hilite_colors = np.concatenate((rgb_vals, alphas), axis=1)
        self.set_color(hilite_colors)
        return self

    def fade(self,depth=0,elev=None,azim=None,ax=None) :
        """
        Reduce surface face opacity based on face position relative
        to the view orientation. 
        
        Parameters
        ----------
        depth : scalar, optional, default: 0
            Minimum opacity to 1 for face opacity from back
            to front face center position.
            Depth value ranges from 0 to 1.

        elev : scalar, optional, default: 30
            Elevation of the axes view.

        azim : scalar, optional, default: -60
            Azimuth of the axes view. 
        
        ax: Matplotlib 3D axes.
            If not None, elev and azim are taken from the ax.
        
        Returns
        -------
        self : surface object

        """

        if ax is None :
            if elev is None: elev = _DFT_VIEW[0]
            if azim is None: azim = _DFT_VIEW[1]
        else :
            elev = ax.elev
            azim = ax.azim
        direction = elev_azim_2vector(elev, azim)
        unitDirection = np.divide( direction, np.linalg.norm(direction) )

        if np.any( (depth<0) or (depth>1) ) :
            raise ValueError('depth values, {}, must be between 0 and 1'.format(depth))

        verts = self.vertexCoor[ np.array(self.fvIndices) ]
        faceCenterCoor = self._get_face_centers(verts)

        self._postProc_surfaceColors()
        colors = self._dfacecolors  # update for v_1.2.0
        if len(faceCenterCoor) == 1 : return self # fade not available for single face line 
        dtprod = np.dot(faceCenterCoor,unitDirection)
        dprange = np.amax(dtprod)-np.amin(dtprod)
        norm_dp = (dtprod - np.amin(dtprod))/dprange
        fadex = (1-depth)*norm_dp + depth

        orig_colors = colors
        fc_less_alpha = orig_colors[:,:3]
        fadex = fadex*orig_colors[:,3]  # update for v_1.2.0
        faded_colors = np.concatenate((fc_less_alpha, fadex[:,np.newaxis] ), axis=1)
        self.set_color(faded_colors) 
        # heuristic approach to hiding face edges.  BANDAID update for v_1.2.0
        alpha = (depth+1)/2  # take the midpoint of the gradient.
        lw =  0.0063*np.exp(4.2*alpha)
        self.set_linewidth(lw)
          
        return self

    def transform(self, rotate=None, scale=1.0, translate=[0,0,0] )  :
        """
        Linear transformation of the surface object.

        Parameters
        ----------
        rotate :  array-like, optional, default: None
            A 3 by 3 rotation matrix.

        scale : 3D array or scalar, optional, default: 1.0
            Multipliers of the object coordinates, in
            xyz coordinates. If a single scalar, three
            multipliers are the same value.
        
        translate : array-like, optional, default: [0,0,0]
            A xyz translation vector of object origin. 

        Returns
        -------
        self : surface object

        """
        if rotate is None : rotate = np.identity(3)
        if isinstance(scale,(int,float)) :
            scale =  np.full(3,scale).astype(float)
        #.....................................    
        def isOrth(M) :
            M = np.array(M)
            tol = 0.00001
            prod = np.dot(M,M.T)
            np.fill_diagonal(prod,0.0)
            return not (np.abs(prod)>tol).any()
        #.....................................    
        trans = lambda V,S,T : np.add(  np.dot( np.multiply(V,S),rotate) ,T)
        #.....................................    
        if not isOrth(rotate) :
            warnings.warn('Transform rotate matrix is NOT Orthoginal.')

        v = trans( self.vertexCoor, scale, translate  )
        if self.phongNorms is not None :    # update for v_1.3.0
            pV = trans( self.phongNorms, scale, translate  )
            self.phongNorms = pV
        if self.vector_field is not None : 
            nScal = np.divide(scale, np.linalg.norm(scale)   )
            self.vector_field = trans( self.vector_field, nScal, [0,0,0]  )
        self.vertexCoor = v
        self._set_geometric_bounds()
        self.set_verts( v[ np.array(self.fvIndices) ] )
        self._surfaceShapeChanged = True
        self.scale = scale
        self.translation = translate
        self.rotation = rotate
        return self

    def get_transformAxis(self, **kargs) :
        """
        Line3DCollection of the 'last' transformed surface coordinate axis.
        
        Parameters
        ----------
        lenmult : scalar or 3D array, optional, default: 1
            Scalar multiplier of the three coordinate axis.

        width : scalar, optional, default: 1.5
            Line width of the coordinate axis.   

        color : a Color or 3D array of Colors, optional, default: ['r','g','b']

        negaxis : boolean, optional, default: False
            If True, include the negative axis, otherwise axis start at origin.

        Returns
        -------
        Line3DCollection object

        """
        kargDft = { 'lenmult':1.0, 'width':1.5, 'color':['r','g','b'], 'negaxis':False }
        KW = KWprocessor(kargDft,'get_transformAxis',**kargs)
        lenmult = KW.getVal('lenmult')
        width = KW.getVal('width')
        color = KW.getVal('color')
        negaxis = KW.getVal('negaxis')
        
        selfProp = {
            'scale' : self.scale,
            'rotation' : self.rotation,
            'translation' : self.translation
        }
        return _transformAxis(selfProp,lenmult, width, color, negaxis)

    def map_geom_from_svd(self, data, pc=None ) :
        """
        Transform surface geometry based on a PCA analysis of a data set.
        Result values contained in the surface property svd_dict.
        
        Parameters
        ----------
        data : N x 3 array of x,y,z data coordinates.
        
        pc : scalar, optional, default: None
            Sets the minimum size of the surface based on the
            relative number of data.  Value range from 0 to 1, 
            with 1 for 100%.  If None, one standard deviation
            of the data sets the surface size.
            
        Returns
        -------
        self : surface object

        """

        if self._isCompositeSurface :
            warnings.warn('Data operation not available for combined shapes.')
            return self
        data_mean = data.mean(axis=0)
        data_centered = data - data_mean

        Ua, sa, Vta = np.linalg.svd(data_centered)
        dataAve = np.mean(sa)       
        scale = sa/dataAve
        tranData = np.divide(np.dot(data_centered,Vta.T),scale )       
        disArr=np.linalg.norm(tranData,axis=1)      
        disArrSum = np.linalg.norm(disArr)
        lendisArr = len(disArr)
        refLen = disArrSum/np.sqrt(lendisArr)        
        s=refLen
        disArr  = np.divide(disArr,refLen)
        
        if pc is not None:
            dlen = len(disArr)
            dindex = int(dlen*pc)
            if dindex >= dlen : dindex = dlen - 1
            if dindex <  0    : dindex = 0
            s = np.sort(disArr)[dindex]
            disArr = disArr/s
            s = s*refLen

        self.transform(Vta,s*scale,data_mean) 
        self._svd_dict = {
            'disarr': disArr, 'sigma': refLen, 'trans': [Vta,s*scale,data_mean]
        }
        self.phongNorms = None  # update for v_1.3.0
        return self

    def contourLines(self,*dist,**kargs) :
        """
        ColorLine3DCollection contour lines on the surface.

        Parameters
        ----------
        dist : floats
            Distances from the origin to the surface of intersection,
            along the direction vector.
              
        direction : array-like, optional, default: [0,0,1]
            A xyz vector normal to the intersection planes
            for planar contours or axial direction of
            cylinders for cylindrical contours.

        name : string identifier (default None)

        color : color (default is the surface colors)

        coor : integer or string indicating the type of contour surface:
            0, p, P, xyz,planar                 - planar (default, planar & polar)
            1, c, C, cylinder,pplar,cylindrical - cylinder (default, cylindrical)
            2, s, S, sphere,spherical           - sphere (default, spherical)
            3, x, X                             - y-z plane
            4, y, Y                             - x-z plane
            5, z, Z                             - x-y plane

        Returns
        -------
        self : line object

        """
        extDft = { 'direction':[0,0,1.0], 'name': None, 'color': None }
        kargDft = { **_COOR_KWARGS, **extDft }
        # reset default 'coor' for cylindrical and spherical coor surfaces.
        coorVal = 0
        if self.coorType == _COORSYS["CYLINDRICAL"] : coorVal = 1
        if self.coorType == _COORSYS["SPHERICAL"] :   coorVal = 2
        if coorVal != 0 :
            temp = copy.deepcopy(_COOR_KWARGS)
            temp['coor'][0] = [coorVal]
            kargDft = { **temp, **extDft }

        KW = KWprocessor(kargDft,'contourLines',**kargs)
        cIndex = KW.getVal('coor')
        direction = KW.getVal('direction')
        name = KW.getVal('name')
        color = KW.getVal('color')

        if direction is None : direction = extDft['direction']
        if cIndex==3 : cIndex, direction = 0, [ 1,0,0 ]
        if cIndex==4 : cIndex, direction = 0, [ 0,1,0 ]
        if cIndex==5 : cIndex, direction = 0, [ 0,0,1 ]

        dotprod = _coor_dotprod(direction,cIndex,self.vertexCoor)
        conLines =  self._contourLineCollection( dotprod,direction,*dist )
        if conLines is None :
            raise ValueError('No surface.contourLines were found.')
        if cIndex == 0 : conLines._planeNormal = direction
        if name is not None  : conLines.name = kargs['name']
        if color is not None : conLines.set_color(kargs['color'])
        return conLines

    def contourLineSet(self,numb=2,**kargs) :
        """
        ColorLine3DCollection of evenly spaced contour lines on the surface.

        Parameters
        ----------
        numb : integer
            Number of contour lines intersecting a surface.
            Numb value must be getter than zero.     
        
        direction : array-like, optional, default: [0,0,1]
            A xyz vector normal to the intersection planes
            for planar contours or axial direction of
            cylinders for cylindrical contours.

        name  : string identifier (default None)

        color : color (default is the surface colors)

        coor : integer or string indicating the type of contour surface:
            0, p, P, xyz,planar                 - planar (default, planar & polar)
            1, c, C, cylinder,pplar,cylindrical - cylinder (default, cylindrical)
            2, s, S, sphere,spherical           - sphere (default, spherical)
            3, x, X                             - y-z plane
            4, y, Y                             - x-z plane
            5, z, Z                             - x-y plane

        Returns
        -------
        self : line object

        """
        extDft = { 'direction':[0,0,1.0], 'name': None, 'color': None }
        kargDft = { **_COOR_KWARGS, **extDft }
        # reset default 'coor' for cylindrical and spherical coor surfaces.
        coorVal = 0
        if self.coorType == _COORSYS["CYLINDRICAL"] : coorVal = 1
        if self.coorType == _COORSYS["SPHERICAL"] :   coorVal = 2
        if coorVal != 0 :
            temp = copy.deepcopy(_COOR_KWARGS)
            temp['coor'][0] = [coorVal]
            kargDft = { **temp, **extDft }
        
        KW = KWprocessor(kargDft,'contourLineSet',**kargs)
        cIndex = KW.getVal('coor')
        direction = KW.getVal('direction')
        name = KW.getVal('name')
        color = KW.getVal('color')

        dotprod = _coor_dotprod(direction,cIndex,self.vertexCoor)  
        maxd = np.amax(dotprod)
        mind = np.amin(dotprod)
        convals = np.linspace(mind,maxd,numb+2)[1:numb+1]
        convals = tuple(convals) 
        conLines = self._contourLineCollection( dotprod,direction,convals )
        if cIndex == 0 : conLines._planeNormal = direction
        if name is not None  : conLines.name = kargs['name']
        if color is not None : conLines.set_color(kargs['color'])
        
        return conLines

    def evert(self) :
        """
        Reverse direction of face normals.
        """
        self.fvIndices = np.flip(self.fvIndices, axis=1)
        if self.phongNorms is not None :    # update for v_1.3.0 
            self.phongNorms = -1.0 * self.phongNorms
        return self

    def domain(self, xlim=None, ylim=None, zlim=None) :
        """
        Set the domain of the surface.

        Used for setting the domain of the base.

        Parameters
        ----------
        xlim, ylim, zlim : 2d arrays
            Minimum and maximum values of the arrays are used to scale
            and translate the surface from an intial domains of
            [ -1, 1 ]
            If set to None, the domains are unchanged.

        Returns
        -------
        self : surface object

        """
        def setDomain(xyz, mx,bx, my,by, mz,bz) :
            x,y,z = xyz
            X = mx*x + bx
            Y = my*y + by
            Z = mz*z + bz
            return X,Y,Z       
        def getCont( lims ) :
            m, b = 1.0, 0.0
            if lims is not None :
                # allow single bndry to be input as a float
                if isinstance(lims, (int,float)) : lims = [-lims,lims]
                m, b = 0.5*(lims[1]-lims[0]) , 0.5*(lims[1]+lims[0])
            return m, b
        mx, bx = getCont(xlim)
        my, by = getCont(ylim)
        mz, bz = getCont(zlim)
        self._map_geom_from_op(lambda A : setDomain(A, mx,bx, my,by, mz,bz),False,None,False)
        return self


class CubicSurface(Surface3DCollection) :
    """
    Cubic 3D surface in Cartesian coordinates with rectangular faces.

    Methods are inherited from the Surface3DCollection.

    """
    
    @staticmethod 
    def _get_cubic_surface_dictionary() :
        """ Constructs the dictionary of base planar surface geometries. """

        surfacesDictionary = {}
        baseVert = []

        baseVert = [ [ 1,-1, -1], [ 1, 1, -1], [-1, 1, -1], [-1,-1, -1],
                     [ 1,-1,  1], [ 1, 1,  1], [-1, 1,  1], [-1,-1,  1] ]

        cubeFaces=(
            [0,1,5,4], [1,2,6,5], [2,3,7,6], [3,0,4,7],
            [4,5,6,7], [0,3,2,1]
        )

        Fc = len(cubeFaces)

        cube =  { 'baseFaceIndices' : cubeFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fc,0,2] }
        surfacesDictionary['cubic'] = cube

        return surfacesDictionary

    _base_surfaces = _get_cubic_surface_dictionary.__func__()
    _default_base = 'cubic'  # default assignment is indicated in __init Docstring.

    @classmethod
    def _Xfev(cls,rez, basetype=None) :
        """
        Calculates the number of faces, edges, and vertices of a CubicSurface.

        Parameters
        ----------
        rez : integer
            Number of recursive subdivisions of a triangulated base faces.
        
        basetype : None
            Placeholder argument has no effect.
        
        Returns
        -------
        list: the number of faces, edges, and vertices.

        """
        F = 6*(4**rez) 
        return [  F, 2*F, F+2  ]

    def __init__(self, rez=0, name=None, **kwargs):
        """
        Create a cubic surface of 2 x 2 x 2 units.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the rectangulated base faces.
            Rez values range from 0 to 7.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.
            
        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.

        """
        # ................................................................
        def midVFun(vectA, vectB) :
            """Mid-point between two points in the xy plane."""
            mid = np.add(vectA,vectB)
            mid = np.multiply(0.5,mid)
            return mid
        # ................................................................

        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        basetype = self._default_base
        baseSurfObj = self._base_surfaces[basetype]
        baseVcoor = baseSurfObj['baseVerticesCoor']
        baseFaceVertexIndices =baseSurfObj['baseFaceIndices']
        
        indexObj,vertexCoor = self._rectangulateBase(rez,baseVcoor,baseFaceVertexIndices, midVFun)
        faceIndices = indexObj['face']
        
        super().__init__(vertexCoor, faceIndices, name, **kwargs)

        self.coorType = _COORSYS["XYZ"]
        self._normal_scale = _normLength( len(faceIndices), 0.61, 1.50 )
        self._rez = rez
        self._basetype = basetype
        return

    def domain(self, xlim, ylim=None, zlim=None) :
        """
        Set the domain of the cubic surface.

        Used for setting the min,max values of the base cube.
        
        Note: this CubicSurface method uses an alternative list of arguments
        from the base class domain method.  Using the xlim arg abbreviates setting
        identical domains for all three axes. 

        Parameters
        ----------
        xlim : a number, list or array, default: 1
            The domain of the function evaluation.  For a number, n,
            the x,y,z axes domains will be [-n,n]. For a 1-dimensional 2-element list,
            [a,b] will be assigned the domain for all 3 axes.  Using a list of list (array),
            as [ [a,b],[c,d],[e,f] ], the domain is assigned individually for each of the
            three coordinate axes.
            Must be a 2D array if ylim or zlim is set.

        ylim, zlim : 2D array, default: None
            Set the ylim and/or zlim.

        Returns
        -------
        self : CubicSurface object

        """
        # the following allows alternative 'arg list' used in super.
        test = ylim is None and zlim is None
        if not test:
            super().domain(xlim,ylim,zlim)
            return self
        domain = _interpret_domain(xlim)    # update for v_1.3.0

        scl = lambda min,max : (max-min)/2
        trn = lambda min,max : (max+min)/2
        scale =     [ scl(d[0],d[1]) for d in domain ]
        translate = [ trn(d[0],d[1]) for d in domain ]
        self.transform(None,scale,translate)
        return self



class PlanarSurface(Surface3DCollection) :
    """
    Flat square 3D surface in Cartesian coordinates.

    Methods are inherited from the Surface3DCollection.
    Planar surface geometries are created in this subclass.

    """
    
    @staticmethod
    def rand(rez=0,seed=None,name=None,kind=None,**kargs) :
        """
        PlanarSurface with random triangular faces.

        Parameters
        ----------
        rez : number, optional, default: 0
            Number of 'recursive' subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.

        seed : integer, optional, default: None
            An initialization seed for random generation.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        kind : string, optional, default: None
            (reserved, not implemented) 

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        Returns
        -------
        surface : PlanarSurface object.

        """

        if not isinstance(rez, (int,float)) :
            raise ValueError('Incorrect rez type, must be a number')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))

        if seed is not None :
            if not isinstance(seed, int) :
                raise ValueError('Incorrect seed type, must ba of type int') 
            np.random.seed(seed)

        N = int( 8*(4**rez)/2 + 4*(2**rez) + 1 )
        x = 2*np.random.rand( N ) - 1
        y = 2*np.random.rand( N ) - 1 
        xydata = np.array([x,y]).T
        tess = spatial.Delaunay(xydata)

        if name is None : name = 'random mesh'
        surface = PlanarSurface(name=name,**kargs)  # set to the class default object
        # reset verts, faceIndices and edges (retain class coordinate system).
        verts = np.zeros( (N,3))
        verts[:,:2] = tess.points
        faceIndices = tess.simplices
        
        surface.vertexCoor = verts
        surface.fvIndices = faceIndices
        surface.evIndices = surface._get_edges_from_faces(faceIndices)
        surface.vfIndicesList = surface._vertexCommonFaceIndices()
        surface.efIndicesList = surface._edgeCommonFaceIndices()
        surface.set_verts( verts[faceIndices] )
        surface._set_geometric_bounds()
        surface._rez = str(rez) if isinstance(rez,int) else '{:.1f}'.format(rez)
        surface._basetype = 'rand_'+str(seed)
        return surface

    @staticmethod
    def _get_planar_surfaces_dictionary() :
        """ Constructs the dictionary of base planar surface geometries. """

        surfacesDictionary = {}

        baseVert = [ [0,0,0],
                     [1,1,0], [-1,1,0], [-1,-1,0], [1,-1,0],
                     [1,0,0], [0,1,0], [-1,0,0], [0,-1,0] ]

        qVert = [ baseVert[i] for i in range(0,5) ]

        quadFaces=( 
            [0,1,2],[0,2,3],[0,3,4],[0,4,1] )
        octFaces1=(
            [0,5,6],[0,6,7],[0,7,8],[0,8,5],
            [5,1,6],[6,2,7],[7,3,8],[8,4,5] )
        octFaces2=(
            [0,5,1],[0,1,6],[0,6,2],[0,2,7],
            [0,7,3],[0,3,8],[0,8,4],[0,4,5] )
        gridFaces=(
            [0,5,1,6], [0,6,2,7], [0,7,3,8], [0,8,4,5]
        )

        Fq = len(quadFaces)
        F1 = len(octFaces1)
        F2 = len(octFaces2)
        Fc = len(gridFaces)
        
        quad =  { 'baseFaceIndices' : quadFaces , 'baseVerticesCoor' : qVert, 'fevParam' : [Fq,2,1] }
        oct1 =  { 'baseFaceIndices' : octFaces1 , 'baseVerticesCoor' : baseVert, 'fevParam' : [F1,4,1] }
        oct2 =  { 'baseFaceIndices' : octFaces2 , 'baseVerticesCoor' : baseVert, 'fevParam' : [F2,4,1] }
        grid =  { 'baseFaceIndices' : gridFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fc,0,2] }

        surfacesDictionary['quad'] = quad
        surfacesDictionary['oct1'] = oct1
        surfacesDictionary['oct2'] = oct2
        surfacesDictionary['squ'] = grid

        # ----------------------------------------------------------
        mr = 0.01
        MR = mr*np.sqrt(2)
        qVertS=copy.copy( surfacesDictionary['quad']['baseVerticesCoor'] )
        qVertS[0] =   [ MR,  0, 0]
        qVertS.append([  0, MR, 0])
        qVertS.append([-MR,  0, 0])
        qVertS.append([0  ,-MR, 0])

        baseVertS = copy.copy(surfacesDictionary['oct1']['baseVerticesCoor'])
        baseVertS[0] =   [ mr, mr, 0]
        baseVertS.append([-mr, mr, 0])
        baseVertS.append([-mr,-mr, 0])
        baseVertS.append([ mr,-mr, 0])

        quadFacesS=( 
            [5,1,2],[6,2,3],[7,3,4],[0,4,1],
            [1,5,0],[2,6,5],[3,7,6],[4,0,7] )

        octFaces1S=(
            [ 0, 5, 6],[ 9, 6, 7],[10, 7, 8],[11, 8, 5],
            [ 5, 1, 6],[ 6, 2, 7],[ 7, 3, 8],[ 8, 4, 5],
            [ 5, 0,11],[ 6, 9, 0],[ 7,10, 9],[ 8,11,10] )

        octFaces2S=(
            [ 0, 5, 1],[ 0, 1, 6],[ 9, 6, 2],[ 9, 2, 7],
            [10, 7, 3],[10, 3, 8],[11, 8, 4],[11, 4, 5],
            [ 5, 0,11],[ 6, 9, 0],[ 7,10, 9],[ 8,11,10] )


        FqS = len(quadFacesS)
        F1S = len(octFaces1S)
        F2S = len(octFaces2S)
        quadS =  { 'baseFaceIndices' : quadFacesS , 'baseVerticesCoor' : qVertS, 'fevParam' : [FqS,2,1] }
        oct1S =  { 'baseFaceIndices' : octFaces1S , 'baseVerticesCoor' : baseVertS, 'fevParam' : [F1S,4,1] }
        oct2S =  { 'baseFaceIndices' : octFaces2S , 'baseVerticesCoor' : baseVertS, 'fevParam' : [F2S,4,1] }

        surfacesDictionary['quad_c'] = quadS
        surfacesDictionary['oct1_c'] = oct1S
        surfacesDictionary['oct2_c'] = oct2S

        return surfacesDictionary
    
    @staticmethod
    def _midVectorFun(vectA, vectB) :
        """Mid-point between two points in the xy plane."""

        mid = np.add(vectA,vectB)
        mid = np.multiply(0.5,mid)
        return mid
    
    _base_surfaces = _get_planar_surfaces_dictionary.__func__()
    _default_base = 'quad'  # default assignment is indicated in __init Docstring.

    @classmethod
    def grid(cls,xdiv, ydiv, basetype='r', minrad=None, name=None, **kwargs) :
        """
        PlanarSurface with X and Y rectangular faces.

        Parameters
        ----------
        xdiv : integer, required.
            Number of X divisions.

        ydiv : integer, required.
            Number of Y divisions.

        basetype : {'d','r'}, optional, default: 'r'
            Starting surface geometries from which the surface is
            constructed indicating face shape.

        minrad : placeholder, not used.

        name : string, optional, default: None.
            Descriptive identifier.

        Raises
        ------
        ValueError
            If xdiv is not an integer greater or equal to 1.
            If ydiv is not an integer greater or equal to 1.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        Returns
        -------
        surface : PlanarSurface object.

        """
        splitSize = None
        minrad = None
        # note: reversal of nang, nrad argsuments compared to other inherited class grid methods.
        nrad, nang = ydiv, xdiv
        surface = cls._square_net(nrad, nang, basetype,minrad, splitSize, name, **kwargs)
        surface._postProc_surfaceColors(True)
        surface._rez = surface._calc_rez( 4 )
        return surface

    @classmethod
    def pntsurf(cls, xyz_points, name=None, **kargs) :
        """
        PlanarSurface from points in Cartesian coordinates.

        Parameters
        ----------
        xyz_points : array, required.
            An N X 3, or N X 4, array of N Cartesian coordinate
            points. For an N X 4 array, the 4th is the vertex
            scalar value.

        name : string, optional, default: None.
            Descriptive identifier for the point surface.

        Returns
        -------
        surface : PlanarSurface object

        """
        dataPts = np.array(xyz_points)
        # update for v_1.3.0
        dataPts, vVal = Surface3DCollection._extract_vals_from_inputCoor(dataPts)

        xyzPts = dataPts.T

        fvIndices = tri.Triangulation(xyzPts[0],xyzPts[1]).triangles
        base_surface = Surface3DCollection(dataPts, fvIndices, color='tan').shade()
        faceIndices = base_surface.fvIndices
        verts = dataPts

        if name is None : name = 'point surface'
        surface = PlanarSurface(name=name,**kargs) 
        surface.vertexCoor = verts
        surface.fvIndices = faceIndices
        surface.evIndices = surface._get_edges_from_faces(faceIndices)
        surface.vfIndicesList = surface._vertexCommonFaceIndices()
        surface.efIndicesList = surface._edgeCommonFaceIndices()
        surface.set_verts( verts[faceIndices] )
        surface._set_geometric_bounds()
        surface._rez = surface._calc_rez(4)
        surface._basetype = 'pntsurf'
        if vVal is not None : surface.set_vertvals(vVal)   # update for v_1.3.0
        return surface

    @classmethod
    def _Xfev(cls,rez, basetype=None) :
        """
        Calculates the number of faces, edges, and vertices of a PlanarSurface.

        Parameters
        ----------
        rez : integer
            Number of recursive subdivisions of a triangulated base faces.
        
        basetype : string
            Starting surface geometries from which the surface would be
            constructed using recursive subdivisions of the triangular faces.
        
        Returns
        -------
        list: the number of faces, edges, and vertices.

        """

        return super()._fev(rez,basetype,cls)

    @staticmethod
    def meshgrid(X,Y,Z,revnormal=False,grid=False, name=None) :
        """
        PlanarSurface object based on a meshgrid.
        
        Parameters
        ----------
        X,Y,Z : N x M arrays
            Cartesian x,y,z coordinate values.

        revnormal : { True, False }, boolean, default: False
            Reverse the face normal directions if True.

        grid : { True, False }, boolean, default: False
            If True, polygon faces are rectangular.  If
            False, polygon faces are triangular.

        name : string, optional, default: None.
            Descriptive identifier for the mesh surface.

        Returns
        -------
        surface : PlanarSurface object
        
        """

        vertCoor = np.reshape(np.stack((X,Y,Z), axis=2 ), (-1,3) )
        lenY, lenX = Z.shape
        mtr = np.reshape( np.arange(lenX*lenY)  , (lenY,lenX))
        A = mtr[:,0:lenX-1][0:lenY-1,:].flatten()
        B = A + 1
        C = A + lenX +1
        D = A + lenX
        if revnormal :
            t = B
            B = D
            D = t
        faceIndices = np.stack( (A,B,C,D) ).T

        if grid : obj = PlanarSurface( basetype='squ')
        else :    obj = PlanarSurface( )
        obj.vertexCoor = vertCoor
        obj.fvIndices = faceIndices
        obj._normal_scale = _normLength( len(faceIndices), 0, 0, 4 )
        obj._set_geometric_bounds()
        obj._surfaceShapeChanged = True        

        v = obj.vertexCoor
        fvi = obj.fvIndices

        if grid :
            obj.vfIndicesList = obj._vertexCommonFaceIndices()
            edgeIndices = obj._get_edges_from_faces(obj.fvIndices)
            obj.evIndices = np.array(edgeIndices)
            obj.efIndicesList = obj._edgeCommonFaceIndices()

        obj.set_verts(v[fvi])

        if not grid : obj.triangulate()
        if name is None : name = 'meshgrid'
        obj.name = name

        return obj

    def __init__(self, rez=0, basetype=None, name=None, **kwargs):
        """
        Create flat square surface of 2 x 2 units.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.
            
        basetype : {'quad','oct1','oct2','squ'}, optional, default: 'quad'
            Starting surface geometries from which the surface is
            constructed using recursive subdivisions of the faces.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.

        """

        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None : basetype = self._default_base
        if basetype not in self._base_surfaces :
            raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(self._base_surfaces)))
        baseSurfObj = self._base_surfaces[basetype]
        baseVcoor = baseSurfObj['baseVerticesCoor']
        baseFaceVertexIndices =baseSurfObj['baseFaceIndices']
        if basetype == 'squ' :
            indexObj,vertexCoor = self._rectangulateBase(rez,baseVcoor,baseFaceVertexIndices, self._midVectorFun)
        else :    
            indexObj,vertexCoor = self._init_triangulateBase(rez,baseVcoor,baseFaceVertexIndices, self._midVectorFun)
        faceIndices = indexObj['face']

        super().__init__(vertexCoor, faceIndices, name, **kwargs)

        self.coorType = _COORSYS["XYZ"]
        self._normal_scale = _normLength( len(faceIndices), 0, 0, 4 )
        self._rez = rez
        self._basetype = basetype
        return

    def scale_dataframe(self,X,Y,Z) :
        """
        Scaling the surface geometry based on a datagrid.

        Parameters
        ----------
        X, Y, Z : N x M arrays
            Minimum and maximum values of the arrays are used to scale
            and translate the surface from an intial domain of
            [ (-1,1), (-1,1), (0,1) ]

        Returns
        -------
        self : PlanarSurface object

        """

        return _scale_dataframe(self,X,Y,Z)

    def domain(self, xlim=None, ylim=None, zcoor=None) :
        """
        Set the domain of the planar surface.

        Used for setting the limits of the base plane.

        Parameters
        ----------
        xlim, ylim : 2d arrays
            Minimum and maximum values of the arrays are used to scale
            and translate the surface from an intial domains of
            [ -1, 1 ]
            If set to None, the domains are unchanged.

        zcoor: number, default : 0
            Z offset of the plane.
            If set to None, the Z offset is unchanged.

        Returns
        -------
        self : PlanarSurface object


        """
        def setDomain(xyz, mx,bx,my,by,mz) :
            x,y,z = xyz
            X = mx*x + bx
            Y = my*y + by
            Z = z if mz is None else np.full(len(z),mz)
            return X,Y,Z       
        def getCont( lims ) :
            m, b = 1.0, 0.0
            if lims is not None :
                # allow single bndry to be input as a float
                if isinstance(lims, (int,float)) : lims = [-lims,lims]
                m, b = 0.5*(lims[1]-lims[0]) , 0.5*(lims[1]+lims[0])
            return m, b
        mx, bx = getCont(xlim)
        my, by = getCont(ylim)
        self.map_geom_from_op(lambda A : setDomain(A, mx,bx,my,by,zcoor))
        return self

class PolarSurface(Surface3DCollection) :
    """
    Flat disk 3D surface in polar coordinates.

    Methods are inherited from the Surface3DCollection.
    Polar surface geometries are created in this subclass.

    """

    @staticmethod
    def rand(rez=0,seed=None,name=None,kind=None,**kargs) :
        """
        PolarSurface with random triangular faces.

        Parameters
        ----------
        rez : number, optional, default: 0
            Number of 'recursive' subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.

        seed : integer, optional, default: None
            An initialization seed for random generation.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        kind : string, optional, default: None
            (reserved, not implemented) 

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        Returns
        -------
        surface : PolarSurface object.

        """
        #.........................................................
        def random_linDist(size) :
            nLoop,nControl,samples = 0, 10**6, []
            while len(samples)<size and nLoop<nControl:
                x=np.random.random_sample()
                prob=x/2
                assert prob>=0 and prob<=1
                if np.random.random_sample() <=prob: samples += [x]
                nLoop+=1
            return np.array(samples)
        #.........................................................

        if not isinstance(rez, (int,float)) :
            raise ValueError('Incorrect rez type, must be a number.')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be a number, 0 <= rez <= {}'.format(rez,_MAXREZ))

        if seed is not None :
            if not isinstance(seed, int) :
                raise ValueError('Incorrect seed type, must ba of type int') 
            np.random.seed(seed)

        N = int( 6*(4**rez)/2 + 3*(2**rez) + 1 )
        r = random_linDist( N )
        t = 2*np.pi*np.random.rand( N )
        z = np.zeros( N )
        data = PolarSurface.coor_convert([r,t,z],True).T
        rtdata = data[:,:2]
        tess = spatial.Delaunay(rtdata)

        if name is None : name = 'random mesh'
        surface = PolarSurface(name=name,**kargs)  # set to the class default object
        # reset verts, faceIndices and edges (retain class coordinate system).
        verts = np.zeros( (N,3))
        verts[:,:2] = tess.points
        faceIndices = tess.simplices
        
        surface.vertexCoor = verts
        surface.fvIndices = faceIndices
        surface.evIndices = surface._get_edges_from_faces(faceIndices)
        surface.vfIndicesList = surface._vertexCommonFaceIndices()
        surface.efIndicesList = surface._edgeCommonFaceIndices()
        surface.set_verts( verts[faceIndices] )
        surface._set_geometric_bounds()
        surface._rez = str(rez) if isinstance(rez,int) else '{:.1f}'.format(rez)
        surface._basetype = 'rand_'+str(seed)
        return surface

    @staticmethod
    def _get_polar_surfaces_dictionary() :
        """Constructs the dictionary of base polar surface geometries."""

        surfacesDictionary = {}
        soff = 0.0001
        x0 = 0.5
        y0 = np.sqrt(3)/2
    
        baseVert = [ [0,0,0],
                     [1,0,0],  [x0,y0,0],   [-x0,y0,0],
                     [-1,0,0], [-x0,-y0,0], [x0,-y0,0] ]

        qVert = [ [0,0,0],
                  [1,0,0], [0,1,0], [-1,0,0], [0,-1,0] ]

        x0S, y0S = soff*y0, soff*x0

        baseVertS = [ [1,-y0S,0], [1, y0S,0],
                      [x0,y0,0], [-x0,y0,0], [-1,0,0], 
                      [-x0,-y0,0], [x0,-y0,0], [ x0S,  y0S,0],
                      [   0, soff,0], [-x0S,  y0S,0],
                      [-x0S, -y0S,0], [   0,-soff,0], [ x0S, -y0S,0] ]

        qVertS = [ [1,-soff,0],
                   [1, soff,0], [0,1,0], [-1,0,0], [0,-1,0],
                   [ soff, soff,0], [-soff ,soff,0], [-soff,-soff,0], [soff,-soff,0] ]
        
        baseFaces=([0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,6],[0,6,1])
        quadFaces=([0,1,2],[0,2,3],[0,3,4],[0,4,1])

        baseFacesS=( [7,1,2], [8,2,3], [9,3,4],  [10,4,5],  [11,5,6], [12,6,0],
                     [2,8,7], [3,9,8], [4,10,9], [5,11,10], [6,12,11] )

        quadFacesS=( [5,1,2],[6,2,3],[7,3,4],[8,4,0],
                    [2,6,5],[3,7,6],[4,8,7] )

        Fs = len(quadFaces)
        Fh = len(baseFaces)
        FsS = len(quadFacesS)
        FhS = len(baseFacesS)
        square =  { 'baseFaceIndices' : quadFaces , 'baseVerticesCoor' : qVert, 'fevParam' : [Fs,2,1] }
        hexag =  { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fh,3,1] }
        squareS =  { 'baseFaceIndices' : quadFacesS , 'baseVerticesCoor' : qVertS, 'fevParam' : [FsS,4.5,1] }
        hexagS =  { 'baseFaceIndices' : baseFacesS , 'baseVerticesCoor' : baseVertS, 'fevParam' : [FhS,6.5,1] }
     
        surfacesDictionary['squ'] = square
        surfacesDictionary['hex'] = hexag
        surfacesDictionary['squ_s'] = squareS
        surfacesDictionary['hex_s'] = hexagS
        # ... 'special case' only for fev method for polar surfaces.
        surfacesDictionary['squ_c'] = { 'fevParam' : [16,5,1] }
        surfacesDictionary['hex_c'] = { 'fevParam' : [12,4,1] }

        return surfacesDictionary

    @staticmethod
    def _midVectorFun(vectA, vectB) :
        """Mid-point at the mid radial position between two points in a plane."""
        mid = np.add(vectA,vectB)
        mid = np.multiply(0.5,mid)
        a = (np.linalg.norm(vectA)+np.linalg.norm(vectB))/2
        b = np.linalg.norm(mid)
        mid = np.multiply(mid,a/b)
        return mid

    @staticmethod
    def _map3Dto2Dsquare(xyz) :
        """Override map surface to rectagular unit coordinates (0,1)."""
        x,y,_ = xyz
        theta = np.arctan2(y,x)
        theta = np.where( theta<0, 2.0*np.pi + theta ,theta)
        a = theta/(2*np.pi)
        r = np.sqrt(x*x+y*y)
        b = 1-r
        a = np.clip(a,0.0,1.0) 
        b = np.clip(b,0.0,1.0) 
        return [a,b]

    @staticmethod
    def coor_convert(xyz, tocart=False) :
        """
        Tranformation between polar and Cartesian coordinates.
        The angular coordinate is in radians with domain 0 to 2pi.

        Parameters
        ----------
        xyz : 3xN array
            N number of vectors in either polar or cartesian coordinates.

        tocart : { True, False }, default: False
            If True, input is polar and output is cartesian.
            If False, input is cartesian and output is polar.

        Returns
        -------
        list: 3xN array
            Transformed coordinates. 

        """

        if tocart :
            r, theta, z = xyz
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            abc = [x,y,z]
        else:
            x,y,z = xyz
            R = np.sqrt( x*x + y*y )
            theta = np.arctan2(y,x)
            theta = np.where(theta<0,theta+2*np.pi,theta)
            abc = [R,theta,z]
        abc = np.array(abc)
        return abc

    _base_surfaces = _get_polar_surfaces_dictionary.__func__()
    _default_base = 'hex'  # default assignment is indicated in __init Docstring.
    _geomType = _COORSYS['POLAR']

    @classmethod
    def _uvw_to_xyz( cls,uvw, rtz) :  
        """Override rotational transformation at a polar coordinate."""
        return super()._disp_to_xyz( uvw, rtz[1] )

    @classmethod
    def grid(cls,nrad, nang, basetype='d', minrad=None, name=None, **kwargs) :
        """
        PolarSurface with axial and radial faces.

        Parameters
        ----------
        nrad : integer, required.
            Number of radial subdivisions.  Minimum value is 1

        nang : integer, required.
            Number of angular subdivisions. Minimum value is 3.

        basetype : {'d','s','q','w','r','x'}, optional, default: 'd'
            Starting surface geometries from which the surface is
            constructed indicating face shape and if split.

        minrad : scalar, optional, default: 0.01
            The minimum distance from the z-axis of any vertex coordinate.
            Use is dependent on basetype.

        name : string, optional, default: None.
            Descriptive identifier.

        Raises
        ------
        ValueError
            If nrad is not an integer greater or equal to 1.
            If nang is not an integer greater or equal to 3.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        Returns
        -------
        surface : PolarSurface object.

        """
        Nfaces = { 'd': 3, 's': 3, 'q':6, 'w':6, 'r':3, 'x':3  }  # 1,3
        splitSize = None
        surface = cls._square_net(nrad, nang, basetype,minrad, splitSize, name, **kwargs)
        surface._postProc_surfaceColors(True)
        surface._rez = surface._calc_rez( 6 )
        return surface

    @classmethod
    def pntsurf(cls, rtz_points, name=None, **kargs) :
        """
        PolarSurface from points in polar coordinates.

        Parameters
        ----------
        rtz_points : array, required.
            An N X 3, or N X 4, array of N polar coordinate
            points. For an N X 4 array, the 4th is the vertex
            scalar value.

        name : string, optional, default: None.
            Descriptive identifier for the point surface.

        Returns
        -------
        surface : PolarSurface object

        """
        dataPts = np.array(rtz_points)
        # update for v_1.3.0
        dataPts, vVal = Surface3DCollection._extract_vals_from_inputCoor(dataPts)

        xyzPts = PolarSurface.coor_convert(dataPts.T,True)

        fvIndices = tri.Triangulation(xyzPts[0],xyzPts[1]).triangles
        base_surface = Surface3DCollection(dataPts, fvIndices, color='tan').shade()
        faceIndices = base_surface.fvIndices
        verts = xyzPts.T

        if name is None : name = 'point surface'
        surface = PolarSurface(name=name,**kargs) 
        surface.vertexCoor = verts
        surface.fvIndices = faceIndices
        surface.evIndices = surface._get_edges_from_faces(faceIndices)
        surface.vfIndicesList = surface._vertexCommonFaceIndices()
        surface.efIndicesList = surface._edgeCommonFaceIndices()
        surface.set_verts( verts[faceIndices] )
        surface._set_geometric_bounds()
        surface._rez = surface._calc_rez(4)
        if vVal is not None : surface.set_vertvals(vVal)   # update for v_1.3.0
        surface._basetype = 'pntsurf'
        return surface

    @classmethod
    def _Xfev(cls,rez, basetype=None, minrad=None) :
        """
        Calculates number of faces, edges, and vertices of a PolarSurface.

        Parameters
        ----------
        rez : integer
            Number of recursive subdivisions of a triangulated base faces.
        
        basetype : string
            Starting surface geometries from which the surface would be
            constructed using recursive subdivisions of the triangular faces.
        
        Returns
        -------
        list: the number of faces, edges, and vertices.
        
        """
        if isinstance(rez,(list,tuple)) : # check for disk...
            nrad = rez[0]
            nang = rez[1]
            # only check is ending in s
            isSplit = basetype[ len(basetype)-1] == 's'
            if minrad is None :
                f = 2*nrad*nang - nang
                e = 3*nrad*nang - nang
                v = nrad*nang + 1
                if isSplit :
                    e += nrad
                    v += nrad
            else :  # any value will trigger this condition.
                f = 2*nrad*nang
                e = 3*nrad*nang + nang
                v = (nrad+1)*nang
                if isSplit :
                    e += nrad    
                    v += nrad + 1    
            return [f,e,v]

        return super()._fev(rez,basetype,cls)

    @staticmethod
    def _flat_split_cyclinder(rez,basetype,minrad, **kwargs) :
        """
        Transform cylindrical base surface to polar coordinates. 

        note :  This method provides an alternative to the 
                poler basetypes of 'hex_s' and 'squ_s'.
                The grids constructed here may provide smoother
                geometries depending on the mapping functions.
        """

        # .................................................
        def cyl2pol(rtz):
            r,t,z = rtz
            b = (1.0+minrad)/2.0
            m = -(1.0-minrad)/2.0
            R = m*z + b
            Z = 0*z
            return R,t,Z
        # .................................................
        if basetype == 'squ_c' : cylbase = 'squ_s' 
        if basetype == 'hex_c' : cylbase = 'tri_s'
        cyl_obj = CylindricalSurface(rez,cylbase, **kwargs).map_geom_from_op(cyl2pol)
        baseFaceVertexIndices = cyl_obj._base_surfaces[cylbase] ['baseFaceIndices']    
        fvIndices = cyl_obj.fvIndices
        vertexCoor = cyl_obj.vertexCoor
        evIndices = cyl_obj.evIndices
        return vertexCoor, fvIndices, evIndices, baseFaceVertexIndices,

    def __init__(self, rez=0, basetype=None, minrad=None, name=None, **kwargs):
        """
        Create flat disk surface of unit radius.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.
            
        basetype : {'squ','hex','squ_s','hex_s','squ_c','hex_c'}, optional, default: 'hex'
            Starting surface geometries from which the surface is
            constructed using recursive subdivisions of the triangular faces.

        minrad : scalar, optional, default: 0.01
            For basetypes 'squ_c' and 'hex_c, the minimum
            distance from the z-axis of any vertex coordinate.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.

        """

        if minrad is None : minrad = _MINRAD
        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None : basetype = self._default_base
        # ... check for 'special case' (ugly code but needed)
        if basetype == 'squ_c' or basetype == 'hex_c' :
            vertexCoor, faceIndices, edgeIndices, bfvi = self._flat_split_cyclinder(rez,basetype,minrad, **kwargs)
            baseFaceVertexIndices = bfvi
        else :
            if basetype not in self._base_surfaces :
                raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(self._base_surfaces)))
            baseSurfObj = self._base_surfaces[basetype]
            baseVcoor = baseSurfObj['baseVerticesCoor']
            baseFaceVertexIndices =baseSurfObj['baseFaceIndices']
            indexObj,vertexCoor = self._init_triangulateBase(rez,baseVcoor,baseFaceVertexIndices, self._midVectorFun)
            faceIndices = indexObj['face']
            edgeIndices = indexObj['edge']

        super().__init__(vertexCoor, faceIndices, name, **kwargs)

        self.coorType = _COORSYS["POLAR"]
        self._normal_scale = _normLength( len(faceIndices), 1.2, 1.55, np.pi )
        self._rez = rez
        self._basetype = basetype
        return

    def domain(self, radius=None, zcoor=None) :
        """
        Set the domain of the polar surface.

        Used for setting the limits of the base polar plane.

        Parameters
        ----------
         radius: number, default : 1
            Radius of the polar surface.
            If set to None, the radius is unchanged.

        zcoor: number, default : 0
            Z offset of the polar plane.
            If set to None, the Z offset is unchanged.
            
        Returns
        -------
        self : PolarSurface object

        """
        def setDomain(rtz, mr, mz) :
            r,t,z = rtz
            R = mr*r
            Z = z if mz is None else np.full(len(z),mz)
            return R,t,Z       
        def getCont(lims) :
            m = 1.0
            if lims is not None : m = lims
            return m
        mr = getCont(radius)
        self.map_geom_from_op(lambda A : setDomain(A, mr,zcoor ))
        return self

class CylindricalSurface(Surface3DCollection) :
    """
    Cylindrical 3D surface in cylindrical coordinates.

    Methods are inherited from the Surface3DCollection.
    Cylindrical surface geometries are created in this subclass.

    """

    @staticmethod
    def rand(rez=0,seed=None,name=None,kind=None,**kargs) :
        """
        CylindricalSurface with random triangular faces.

        Parameters
        ----------
        rez : number, optional, default: 0
            Number of 'recursive' subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.

        seed : integer, optional, default: None
            An initialization seed for random generation.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        kind : string, optional, default: None
            (reserved, not implemented) 

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        Returns
        -------
        surface : CylindricalSurface object.

        """
        #.........................................................
        def random_sinDist(size) :
            nLoop,nControl,samples = 0, 10**6, []
            while len(samples)<size and nLoop<nControl:
                x=np.random.random_sample()
                prob=np.sin(x*np.pi)/2
                assert prob>=0 and prob<=1
                if np.random.random_sample() <=prob: samples += [x]
                nLoop+=1
            return np.array(samples)
        #.........................................................

        if not isinstance(rez, (int, float)) :
            raise ValueError('Incorrect rez type, must be a number.')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be a number, 0 <= rez <= {}'.format(rez,_MAXREZ))

        if seed is not None :
            if not isinstance(seed, int) :
                raise ValueError('Incorrect seed type, must ba of type int') 
            np.random.seed(seed)

        N = int(12*(4**rez)/2 + 3*(2**rez))
        r = np.ones( N )
        t = np.random.uniform(0,2*np.pi, N)
        z = np.random.uniform(-1,1, N)
        dataPts = np.array( [r,t,z] ).T
        if name is None : name = 'random mesh'
        surface = CylindricalSurface.pntsurf(dataPts, name=None, **kargs)
        surface._rez = str(rez) if isinstance(rez,int) else '{:.1f}'.format(rez)
        surface._basetype = 'rand_'+str(seed)
        return surface

    @staticmethod
    def _get_cylindrical_surfaces_dictionary() :
        """Constructs the dictionary of base cylindrical surface geometries. """

        surfacesDictionary = {}
        def construct_tricylinder() :
            y0 = np.sqrt(3.0)/2.0
            x0 = 0.5

            baseVert = [ [1,0,0],   [-x0,-y0,0], [-x0,y0,0],
                         [-1,0,1],  [x0,-y0,1],  [x0,y0,1],
                         [-1,0,-1], [x0,-y0,-1], [x0,y0,-1] ]

            baseFaces=( 
                [0,2,5],[2,1,3],[1,0,4],
                [0,5,4],[2,3,5],[1,4,3], 
                [0,8,2],[2,6,1],[1,7,0],
                [0,7,8],[2,8,6],[1,6,7])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,3,0] }

        def construct_tricylinder_v() :
            vfi_d = construct_tricylinder()
            baseVert =  vfi_d['baseVerticesCoor']
            baseFaces = vfi_d['baseFaceIndices']
            baseVert = baseVert + [[0,0,1],[0,0,-1]]
            baseFaces = list(baseFaces)
            addFaces = [ 
                [3,4,9],[4,5,9],[5,3,9],
                [6,8,10],[8,7,10],[7,6,10] ]
            baseFaces = baseFaces + addFaces
            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
            
        def construct_tri2cylinder() :
            y0 = np.sqrt(3.0)/2.0
            x0 = 0.5 

            baseVert = [ [-1,0,0], [x0,-y0,0],   [x0,y0,0],
                         [1,0,1],  [-x0,-y0,1],  [-x0,y0,1],
                         [1,0,-1], [-x0,-y0,-1], [-x0,y0,-1] ]

            baseFaces=( 
                [0,4,5],[1,3,4],[2,5,3],
                [0,8,7],[1,7,6],[2,6,8],
                [0,7,4],[1,6,3],[2,8,5],
                [0,5,8],[1,4,7],[2,3,6])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,3,0] }

        def construct_tri2cylinder_v() :
            vfi_d = construct_tri2cylinder()
            baseVert =  vfi_d['baseVerticesCoor']
            baseFaces = vfi_d['baseFaceIndices']
            baseVert = baseVert + [[0,0,1],[0,0,-1]]
            baseFaces = list(baseFaces)
            addFaces = [ 
                [3,9,4],[4,9,5],[5,9,3],
                [6,10,8],[7,10,6],[8,10,7] ]
            baseFaces = baseFaces + addFaces
            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }

        def construct_tri2cylinder_sliced() :
            y0 = np.sqrt(3.0)/2.0
            x0 = 0.5   
            yof = 0.0001

            baseVert = [ [-1,0,0],   [x0,-y0,0],   [x0,y0,0],
                         [1,yof,1],  [-x0,-y0,1],  [-x0,y0,1],
                         [1,yof,-1], [-x0,-y0,-1], [-x0,y0,-1],
                         [1,-yof,1], [1,-yof,-1] ]

            baseFaces=( 
                [0,4,5],[1,9,4],[2,5,3],
                [0,8,7],[1,7,10],[2,6,8],
                [0,7,4],[1,10,9],[2,8,5],

                [0,5,8],[1,4,7],[2,3,6])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,4,1] }

        def construct_tri2cylinder_sliced_v() :
            vfi_d = construct_tri2cylinder_sliced()
            baseVert =  vfi_d['baseVerticesCoor']
            baseFaces = vfi_d['baseFaceIndices']
            baseVert = baseVert + [[0,0,1],[0,0,-1]]
            baseFaces = list(baseFaces)
            addFaces = [ 
                [3,5,11],[5,4,11],[4,9,11], 
                [8,6,12],[7,8,12],[10,7,12] ]
            baseFaces = baseFaces + addFaces
            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }

        def construct_quadcylinder() :
            xy0 = 1.0/np.sqrt(2.0)

            baseVert = [ [1,0,0],      [0,-1,0],      [-1,0,0],       [0,1,0],
                         [xy0,xy0,1],  [-xy0,xy0,1],  [-xy0,-xy0,1],  [xy0,-xy0,1],
                         [xy0,xy0,-1], [-xy0,xy0,-1], [-xy0,-xy0,-1], [xy0,-xy0,-1] ]

            baseFaces=( 
                [1,0,7],[2,1,6],[3,2,5],[0,3,4],
                [4,3,5],[5,2,6],[6,1,7],[7,0,4],
                [0,8,3],[3,9,2],[2,10,1],[1,11,0],
                [0,11,8],[3,8,9],[2,9,10],[1,10,11] )

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,4,0] }

        def construct_quadcylinder_v() :
            vfi_d = construct_quadcylinder()
            baseVert =  vfi_d['baseVerticesCoor']
            baseFaces = vfi_d['baseFaceIndices']
            baseVert = baseVert + [[0,0,1],[0,0,-1]]
            baseFaces = list(baseFaces)
            addFaces = [ 
                [4,5,12],[5,6,12],[6,7,12],[7,4,12],
                [9,8,13],[10,9,13],[11,10,13],[8,11,13] ]
            baseFaces = baseFaces + addFaces
            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }

        def construct_quad2cylinder() :
            xy0 = 1.0/np.sqrt(2.0)

            baseVert = [ [xy0,xy0,0], [-xy0,xy0,0], [-xy0,-xy0,0], [xy0,-xy0,0],
                         [1,0,1],     [0,-1,1],     [-1,0,1],      [0,1,1],
                         [1,0,-1],    [0,-1,-1],    [-1,0,-1],     [0,1,-1] ]

            baseFaces=( 
                [0,7,4],[1,6,7],[2,5,6],[3,4,5],
                [0,8,11],[1,11,10],[2,10,9],[3,9,8],
                [0,4,8],[1,7,11],[2,6,10],[3,5,9],
                [0,11,7],[1,10,6],[2,9,5],[3,8,4] )

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,4,0] }

        def construct_quad2cylinder_v() :
            vfi_d = construct_quad2cylinder()
            baseVert =  vfi_d['baseVerticesCoor']
            baseFaces = vfi_d['baseFaceIndices']
            baseVert = baseVert + [[0,0,1],[0,0,-1]]
            baseFaces = list(baseFaces)
            addFaces = [ 
                [5,4,12],[6,5,12],[7,6,12],[4,7,12],
                [8,9,13],[9,10,13],[10,11,13],[11,8,13] ]
            baseFaces = baseFaces + addFaces
            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }

        def construct_quad2cylinder_sliced() :
            xy0 = 1.0/np.sqrt(2.0)
            yof = 0.0001

            baseVert = [ [xy0,xy0,0], [-xy0,xy0,0], [-xy0,-xy0,0], [xy0,-xy0,0],
                         [1,yof,1],   [0,-1,1],     [-1,0,1],       [0,1,1],
                         [1,yof,-1],  [0,-1,-1],    [-1,0,-1],      [0,1,-1],
                         [1,-yof,1],  [1,-yof,-1] ]

            baseFaces=( 
                [0,7,4],[1,6,7],[2,5,6],[3,12,5],
                [0,8,11],[1,11,10],[2,10,9],[3,9,13],
                [0,4,8],[1,7,11],[2,6,10],[3,5,9],
                [0,11,7],[1,10,6],[2,9,5],[3,13,12] )

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,5,1] }

        def construct_quad2cylinder_sliced_v() :
            vfi_d = construct_quad2cylinder_sliced()
            baseVert =  vfi_d['baseVerticesCoor']
            baseFaces = vfi_d['baseFaceIndices']
            baseVert = baseVert + [[0,0,1],[0,0,-1]]
            baseFaces = list(baseFaces)
            addFaces = [ 
                [4,7,14],[7,6,14],[6,5,14],[5,12,14],
                [11,8,15],[10,11,15],[9,10,15],[13,9,15] ]
            baseFaces = baseFaces + addFaces
            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }

        surfacesDictionary['tri'] = construct_tricylinder()
        surfacesDictionary['squ'] = construct_quadcylinder()
        surfacesDictionary['tri2'] = construct_tri2cylinder()
        surfacesDictionary['squ2'] = construct_quad2cylinder()
        surfacesDictionary['tri_s'] = construct_tri2cylinder_sliced()
        surfacesDictionary['squ_s'] = construct_quad2cylinder_sliced()

        surfacesDictionary['tri_v'] = construct_tricylinder_v()
        surfacesDictionary['squ_v'] = construct_quadcylinder_v()
        surfacesDictionary['tri2_v'] = construct_tri2cylinder_v()
        surfacesDictionary['squ2_v'] = construct_quad2cylinder_v()
        surfacesDictionary['tri_sv'] = construct_tri2cylinder_sliced_v()
        surfacesDictionary['squ_sv'] = construct_quad2cylinder_sliced_v()
        return surfacesDictionary
    
    @staticmethod
    def _midVectorFun(vectA, vectB) :  
        """Mid-point between two points on a unit cylinder."""
        mid = np.add(vectA,vectB)
        mid = np.multiply(0.5,mid)
        z = np.dot(mid,[0,0,1])
        v = np.subtract(mid,[0,0,z])
        xy =  v/np.linalg.norm(v)
        # xy is unit vector at cylinder radius, but
        # is average of radius for top & bottom faces.
        zA = np.dot(vectA,[0,0,1])
        vA = np.subtract(vectA,[0,0,zA])
        zB = np.dot(vectB,[0,0,1])
        vB = np.subtract(vectB,[0,0,zB])
        a = ( np.linalg.norm(vA) + np.linalg.norm(vB) )/2.0
        xy =  a*xy
        return np.add(xy,[0,0,z])

    @staticmethod
    def _map3Dto2Dsquare(xyz) :
        """Override map surface to rectagular unit coordinates (0,1)."""
        x,y,z = xyz
        theta = np.arctan2(y,x)
        theta = np.where( theta<0, 2.0*np.pi + theta ,theta)
        a = theta/(2*np.pi)
        b =  z/2 + 0.5 
        a = np.clip(a,0.0,1.0) 
        b = np.clip(b,0.0,1.0) 
        return [a,b]

    @staticmethod
    def coor_convert(xyz, tocart=False) :
        """
        Transformation between cylindrical and Cartesian coordinates.
        The angular coordinate is in radians with domain 0 to 2pi.

        Parameters
        ----------
        xyz : 3xN
            N number of vectors in either cylindrical or cartesian coordinates,
            including the z coordinate( which remains unchanged).

        tocart : { True, False }, default: False
            If True, input is cylindrical and output is cartesian.
            If False, input is cartesian and output is cylindrical.

        Returns
        -------
        list: 3xN array
            Transformed coordinates. 

        """

        if tocart :
            r, theta, z = xyz
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            abc = [x,y,z]
        else:
            x,y,z = xyz
            R = np.sqrt( x*x + y*y )
            theta = np.arctan2(y,x)
            theta = np.where(theta<0,theta+2*np.pi,theta)
            abc = [R,theta,z]
        abc = np.array(abc)
        return abc

    _base_surfaces = _get_cylindrical_surfaces_dictionary.__func__()
    _default_base = 'tri'  # default assignment is indicated in __init Docstring.
    _geomType = _COORSYS['CYLINDRICAL']

    @classmethod
    def _cyl_vol_square_net(cls, numZ, numU, basetype,minRad, splitSize, numR, name=None, **kwargs) :
        """
        method only called by the CylinderSurface classmethod grid ONLY for volume cylinder objects.
        """
        # DevNote:  this method uses the same algorithm as the base class _square_net classmethod.
        #           This was done to simplify the code for this special case rather than
        #           add additional 'convoluted' logic to the _square_net method.  
        #           For this case, the algorithm combines one cylindrical and two polar grids.

        # Note: numX and num_X are the number of divisions and the number of coor positions.
        if numR is None :  numR = int((numZ+1)/2)
        if basetype is None : basetype = 'd'
        if minRad is None : minRad = _MINRAD
        if splitSize is None : splitSize = _SPLITSIZE
        
        basetype = basetype.lower()
        baselist = ['d','s','q','w','r','x']
        if basetype not in baselist :
            raise ValueError("Grid basetype '{}' is not recognized. Possible values are: {}.".format(basetype,list(baselist)))

        isSplit = (basetype == 's') or (basetype == 'w') or (basetype == 'x')
        needsTriag = not ( (basetype == 'r') or (basetype == 'x')  )
        hasCenter = (basetype == 'd') or (basetype == 's')

        # T angle ranges   .....
        num_U = numU+1 if isSplit else numU
        min_U = 0.0
        max_U = (1-1/num_U)*2*np.pi
        if isSplit : max_U = (2-splitSize)*np.pi
        coor_U = np.linspace(min_U,max_U,num_U)

        # R disk ranges ........
        if hasCenter :
            num_R,inn_R = numR-1, 1.0/numR
            out_R = 1.0-(1.0/numR)
        else :
            num_R,inn_R = numR,   minRad
            out_R = 1.0-(1.0/numR) + minRad/numR
        
        coor_R_bot = np.linspace(inn_R,out_R, num_R)
        coor_R_mid = np.full(numZ+1,1.0)
        coor_R_top = np.linspace(out_R,inn_R, num_R)
        num_Z = numZ+1
        coor_Z_bot = np.full(num_R,-1.0)
        coor_Z_mid = np.linspace(-1.0,1.0,num_Z)
        coor_Z_top = np.full(num_R, 1.0)

        coor_R = np.concatenate( (coor_R_bot,coor_R_mid,coor_R_top),axis=0)
        coor_Z = np.concatenate( (coor_Z_bot,coor_Z_mid,coor_Z_top),axis=0)
        
        numV = 2*numR + numZ
        num_V = 2*num_R + num_Z
        totVerts = num_U*num_V

        U_coor = np.tile(coor_U,num_V)
        R_coor = np.tile(coor_R,(num_U,1))
        R_coor = np.ravel(R_coor,order='F')
        Z_coor = np.tile(coor_Z,(num_U,1))
        Z_coor = np.ravel(Z_coor,order='F')

        abc = np.array([R_coor,U_coor,Z_coor])
        xyz = CylindricalSurface.coor_convert(abc,True)
        verts = np.array(xyz).T
        numVt = num_V-1 if hasCenter else numV
        faceIndices = cls._grid_face_indices(numU,numVt,isSplit,needsTriag)

        if hasCenter :
            # add top triangular faces ...
            topVert = np.array([[0.0,0.0,1]])
            topIndex = totVerts
            startIndex = topIndex - num_U
            endIndex = topIndex-2 if isSplit else topIndex-1
            a = np.linspace(startIndex,endIndex,numU,dtype=int)
            b = a+1 if isSplit else np.append( a[1:], [startIndex] )
            c = np.full(numU,topIndex)
            abc = np.array([a,b,c]).T
            verts = np.append(verts,topVert, axis=0)
            faceIndices = np.append(faceIndices,abc, axis=0)

            # add bottom triangular facess ...
            btmVert = np.array([[0.0,0.0,-1]])
            btmIndex = totVerts+1
            startIndex,endIndex = 0,numU-1
            a = np.linspace(startIndex,endIndex,numU,dtype=int)
            b = a+1 if isSplit else np.append( a[1:], [startIndex] )
            c = np.full(numU,btmIndex)
            bac = np.array([b,a,c]).T  # bac instead for abc for outward normals.
            verts = np.append(verts,btmVert, axis=0)
            faceIndices = np.append(faceIndices,bac, axis=0)

        surface = cls(**kwargs)  # set to the class default object
        # reset verts, faceIndices and edges (retain class coordinate system).
        surface.vertexCoor = verts
        surface.fvIndices = faceIndices
        surface.evIndices = surface._get_edges_from_faces(faceIndices)
        surface.vfIndicesList = surface._vertexCommonFaceIndices()
        surface.efIndicesList = surface._edgeCommonFaceIndices()
        surface.set_verts( verts[faceIndices] )
        surface._set_geometric_bounds()
        surface.name = name
        surface._basetype = 'grid_'+basetype+'v'
        return surface

    @classmethod
    def _uvw_to_xyz(cls, uvw, rtz) :
        """Override rotational tranformation at a cylindrical coordinate.""" 
        return super()._disp_to_xyz( uvw, rtz[1], None )

    @classmethod
    def grid(cls, nver ,nang, basetype='d',minrad=None, numR=None, name=None, **kwargs) :
        """
        CylindricalSurface with axial and vertical faces.

        Parameters
        ----------
        nver : integer, required.
            Number of vertical subdivisions.  Minimum value is 1

        nang : integer, required.
            Number of angular subdivisions. Minimum value is 3.
        
        basetype : {'d','q','r','s','w','x', 'dv','qv','rv','sv','wv','xv'}, optional, default: 'd'
            Starting surface geometries from which the surface is
            constructed indicating face shape and if split.

        minrad : scalar, optional, default: 0.01
            The minimum distance from the z-axis of any vertex coordinate.
            Use is dependent on basetype.

        numR : integer, optional, default: None
            Only applicable for volume cylinder geometries.
            Number of radial subdivisions for volume cylinders.
            When None, default is half the number of vertical subdivisions.

        name : string, optional, default: None.
            Descriptive identifier.

        Raises
        ------
        ValueError
            If nver is not an integer greater or equal to 1.
            If nang is not an integer greater or equal to 3.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        Returns
        -------
        surface : CylindricalSurface object.

        """
        Nfaces = { 'd': 6,   's': 6,   'q':6,   'w':6,   'r':3,  'x':3,
                   'dv': 12, 'sv': 12, 'qv':18, 'wv':18, 'rv':9, 'xv':9  }  # 1,3
        splitSize = None  # FutDev: currently, not user accessable.
        dftNfaces = 12

        basetype = basetype.lower()
        if basetype.endswith('v') :
            dftNfaces = 18
            subtype = basetype[0]
            surface = CylindricalSurface._cyl_vol_square_net( nver ,nang, subtype, minrad, splitSize, numR, name, **kwargs)
        else : 
            surface = cls._square_net(nver ,nang,basetype,minrad,splitSize, name, **kwargs)
        surface._postProc_surfaceColors(True)
        surface._rez = surface._calc_rez( dftNfaces )
        return surface       

    @classmethod
    def pntsurf(cls, rtz_points, name=None, **kargs) :
        """
        CylindricalSurface from points in cylindrical coordinates.

        Parameters
        ----------
        rtz_points : array, required.
            An N X 3, or N X 4, array of N cylindrical coordinate
            points. For an N X 4 array, the 4th is the vertex
            scalar value.

        name : string, optional, default: None.
            Descriptive identifier for the point surface.

        Returns
        -------
        surface : CylindricalSurface object

        """
        dataPts = np.array(rtz_points)
        # update for v_1.3.0
        dataPts, vVal = Surface3DCollection._extract_vals_from_inputCoor(dataPts)

        cylPts = np.array(dataPts.T)
        cylPts[0] = np.amax(cylPts[0])
        xyz_cylPts = CylindricalSurface.coor_convert(cylPts,True)
        minZ,maxZ = np.amin(cylPts[2]), np.amax(cylPts[2])
        offsetZ = (maxZ+minZ)/2
        xyz_cylPts[2] = xyz_cylPts[2] - offsetZ

        sphPts = SphericalSurface.coor_convert(xyz_cylPts,False)
        rmax = np.amax(sphPts[0])
        sphPts[0] = rmax
        xyz_sphPts = SphericalSurface.coor_convert(sphPts,True)
        ptTop,ptBtm = [0,0,rmax], [0,0,-rmax]
        xyz_sphPts = np.append(xyz_sphPts.T,[ptTop,ptBtm],axis=0).T

        # note: the following surface method chull corrects face normal
        # orientations so spatial.ConvexHull was not directly used.
        hull_surface = Surface3DCollection.chull(xyz_sphPts.T)

        verts = CylindricalSurface.coor_convert(dataPts.T,True).T
        # remove faces attached to the two top & bottom vertices.
        N = len(dataPts)
        fvIndices = hull_surface.fvIndices
        maxIndices = np.max(fvIndices,axis=1)
        shouldkeep = np.where(maxIndices < N)
        faceIndices = fvIndices[shouldkeep]

        if name is None : name = 'point surface'
        surface = CylindricalSurface(name=name,**kargs) 
        surface.vertexCoor = verts
        surface.fvIndices = faceIndices
        surface.evIndices = surface._get_edges_from_faces(faceIndices)
        surface.vfIndicesList = surface._vertexCommonFaceIndices()
        surface.efIndicesList = surface._edgeCommonFaceIndices()
        surface.set_verts( verts[faceIndices] )
        surface._set_geometric_bounds()
        surface._rez = surface._calc_rez(12)
        surface._basetype = 'pntsurf'
        if vVal is not None : surface.set_vertvals(vVal)   # update for v_1.3.0
        return surface

    @classmethod
    def _Xfev(cls,rez, basetype=None, minrad=None) :
        """
        Calculates number of faces, edges, and vertices of a CylindricalSurface.

        Parameters
        ----------
        rez : integer
            Number of recursive subdivisions of a triangulated base faces.
        
        basetype : string
            Starting surface geometries from which the surface would be
            constructed using recursive subdivisions of the triangular faces.
        
        Returns
        -------
        list: the number of faces, edges, and vertices.

        """
        if isinstance(rez,(list,tuple)) : # check for ring...
            nver = rez[0]
            nang = rez[1]
            # only check is ending in s
            isSplit = basetype[ len(basetype)-1] == 's'
            f = 2*nver*nang
            e = 3*nver*nang + nang
            v = (nver + 1)*nang
            if isSplit :
                e += nver
                v += nver + 1
            return [f,e,v]           

        return super()._fev(rez,basetype,cls)
               
    def __init__(self, rez=0, basetype=None, name=None, **kwargs):
        """
        Create cylindrical surface of unit radius and length 2.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.
            
        basetype : {'tri','squ','tri2','squ2','tri_s','squ_s', 'tri_v','squ_v','tri2_v','squ2_v','tri_sv','squ_sv'}, optional, default: 'tri'
            Starting surface geometries from which the surface is
            constructed using recursive subdivisions of the triangular faces.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        """

        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None : basetype = self._default_base
        if basetype not in self._base_surfaces :
            raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(self._base_surfaces)))
        baseSurfObj = self._base_surfaces[basetype]
        baseVcoor = baseSurfObj['baseVerticesCoor']
        baseFaceVertexIndices =baseSurfObj['baseFaceIndices']
        indexObj,vertexCoor = self._init_triangulateBase(rez,baseVcoor,baseFaceVertexIndices, self._midVectorFun)
        faceIndices = indexObj['face']

        super().__init__(vertexCoor, faceIndices, name, **kwargs)

        self.coorType = _COORSYS["CYLINDRICAL"]
        self._normal_scale =  _normLength( len(faceIndices), 1.14, 1.31)
        self._rez = rez
        self._basetype = basetype
        self.disp_Vector = np.array( [1,1,0] )
        return None

    def domain(self, radius=None, zlim=None) :
        """
        Set the domain of the cylindrical surface.

        Used for setting the limits of the base cylinder.

        Parameters
        ----------
        radius: number, default : 1
            Radius of the cylindrical surface.
            If set to None, the domain is unchanged.

        zlim : 2d array
            Minimum and maximum values of the arrays are used to scale
            and translate the surface from an intial domain of
            [ -1, 1 ]
            If set to None, the Z-domain is unchanged.

        Returns
        -------
        self : CylindricalSurface object

        """
        def setDomain(rtz, mr, mz,bz) :
            r,t,z = rtz
            R = mr*r
            Z = mz*z + bz
            return R,t,Z       
        def getCont_r(lims) :
            m = 1.0
            if lims is not None : m = lims
            return m
        def getCont( lims ) :
            m, b = 1.0, 0.0
            if lims is not None :
                # allow single bndry to be input as a float
                if isinstance(lims, (int,float)) : lims = [-lims,lims]
                m, b = 0.5*(lims[1]-lims[0]) , 0.5*(lims[1]+lims[0])
            return m, b
        mr = getCont_r(radius)
        mz, bz = getCont(zlim)
        self.map_geom_from_op(lambda A : setDomain(A, mr, mz,bz))
        return self

class SphericalSurface(Surface3DCollection) :
    """
    Spherical 3D surface in spherical coordinates.

    Methods are inherited from the Surface3DCollection.
    Spherical surface geometries are created in this subclass.

    """

    @staticmethod
    def rand(rez=0,seed=None,name=None,kind=None,**kargs) :
        """
        SphericalSurface with random triangular faces.

        Parameters
        ----------
        rez : number, optional, default: 0
            Number of 'recursive' subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.

        seed : integer, optional, default: None
            An initialization seed for random generation.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        kind : string, optional, default: None
            (reserved, not implemented) 

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        Returns
        -------
        surface : SphericalSurface object.

        """
        #.........................................................
        def random_sinDist(size) :
            nLoop,nControl,samples = 0, 10**6, []
            while len(samples)<size and nLoop<nControl:
                x=np.random.random_sample()
                prob=np.sin(x*np.pi)/2
                assert prob>=0 and prob<=1
                if np.random.random_sample() <=prob: samples += [x]
                nLoop+=1
            return np.array(samples)
        #.........................................................

        if not isinstance(rez, (int, float)) :
            raise ValueError('Incorrect rez type, must be a number.')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be a number, 0 <= rez <= {}'.format(rez,_MAXREZ))

        if seed is not None :
            if not isinstance(seed, int) :
                raise ValueError('Incorrect seed type, must ba of type int') 
            np.random.seed(seed)

        N = int(20*(4**rez)/2 + 2)
        r = np.ones( N )
        t = 2*np.pi*np.random.rand( N )
        p = np.pi*random_sinDist( N )
        data = SphericalSurface.coor_convert([r,t,p],True).T
        hullsurf = Surface3DCollection.chull(data)

        if name is None : name = 'random mesh'
        surface = SphericalSurface(name=name,**kargs)  # set to the class default object
        # reset verts, faceIndices and edges (retain class coordinate system).
        verts = hullsurf.vertexCoor
        faceIndices = hullsurf.fvIndices

        surface.vertexCoor = verts
        surface.fvIndices = faceIndices
        surface.evIndices = surface._get_edges_from_faces(faceIndices)
        surface.vfIndicesList = surface._vertexCommonFaceIndices()
        surface.efIndicesList = surface._edgeCommonFaceIndices()
        surface.set_verts( verts[faceIndices] )
        surface._set_geometric_bounds()
        surface._rez = str(rez) if isinstance(rez,int) else '{:.1f}'.format(rez)
        surface._basetype = 'rand_'+str(seed)
        return surface

    @staticmethod
    def platonic(rez=0,basetype=None,name=None, **kwargs) :
        """
        Platonic Solid surfaces.
        
        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.
            
        basetype : {'tetra','octa','icosa','cube','dodeca'}, optional, default: 'icosa'
            Starting surface geometries from which the surface is
            constructed using recursive subdivisions of the triangular faces.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Returns
        -------
        SphericalSurface object

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.SphericalSurface'.

        """
        #.....................................
        def fev_cube_a() :
            v = [ [-1,-1,-1],[-1, 1,-1],[ 1, 1,-1],[ 1,-1,-1],
                  [-1,-1, 1],[-1, 1, 1],[ 1, 1, 1],[ 1,-1, 1] ]
            f = [ [0,1,2,3],[4,7,6,5],[0,3,7,4],[3,2,6,7],[2,1,5,6],[1,0,4,5]]
            e = [ [3,2],[2,1],[1,0],[0,3], [7,6],[6,5],[5,4],[4,7], [0,4],[3,7],[2,6],[1,5]]
            return np.array(v),np.array(f),np.array(e)
        def _preTriangulateDodoc(v,f,e,fc) :
            # (special Case)
            # Only used for rez=1 of a platonic dodecahedron. Each pentagon
            # is broken into 5 isoceles triangles from a centered vertex.
            subFaceIndices,subEdgeIndices = [],[]
            newCoors,newColor = [],[]
            nvInx = len(v)
            for i, face in enumerate(f) :
                for n0 in range(5) :
                    n1 = (n0+1)%5
                    subFaceIndices.append( [face[n0],face[n1],nvInx+i] )
                    subEdgeIndices.append( [face[n0],nvInx+i] )
                newCoors.append( np.sum(v[face], axis=0)/5 )
                newColor += [fc[i]]*5
            newverts = np.append(v,newCoors,axis=0) 
            newfvInx = np.array(subFaceIndices) 
            newEvInx = np.append(e,subEdgeIndices,axis=0) 
            newColor = np.array(newColor)
            return newverts,newfvInx,newEvInx,newColor
        #.....................................
        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None : basetype = SphericalSurface._default_base

        baseDict = { 'tetra':"tetrahedron",'octa':"octahedron",
            'icosa':"icosahedron",'cube':"cube",'dodeca':"dodecahedron",
            'cube_a':"cube_align"  }
        baseCheck = baseDict.keys()

        if basetype not in baseCheck :
            raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,baseCheck) )

        if basetype == 'cube_a' :
            surface = SphericalSurface(basetype='cube',name=name, **kwargs)
        else :    
            surface = SphericalSurface(basetype=basetype,name=name, **kwargs)

        surface._rez = rez
        if basetype == 'cube_a' :
            v,f,e = fev_cube_a()
            surface.vertexCoor = v
            surface.fvIndices = f
            surface.evIndices = e
            surface.vfIndicesList = surface._vertexCommonFaceIndices()
            surface.efIndicesList = surface._edgeCommonFaceIndices()
            surface.set_verts( v[f] )
        if basetype == 'cube' :
            v,f,e = SphericalSurface.get_cube()
            surface.vertexCoor = v
            surface.fvIndices = f
            surface.evIndices = e
            surface.vfIndicesList = surface._vertexCommonFaceIndices()
            surface.efIndicesList = surface._edgeCommonFaceIndices()
            surface.set_verts( v[f] )
        if basetype == 'dodeca' :
            v,f,e = SphericalSurface.get_dodecahedron()
            if rez>0 :
                fc = surface._dfacecolors  # update for v_1.2.0
                v,f,e,fc = _preTriangulateDodoc(v,f,e,fc)
                surface.set_facecolor(fc)
                rez -= 1
            surface.vertexCoor = v
            surface.fvIndices = f
            surface.evIndices = e
            surface.vfIndicesList = surface._vertexCommonFaceIndices()
            surface.efIndicesList = surface._edgeCommonFaceIndices()
            surface.set_verts( v[f] )

        surface._postProc_surfaceColors(True)
        if rez > 0 :
            surface.triangulate(rez)
      
        surface._set_geometric_bounds()
        surface._basetype = basetype + '*'   
        surface.name = baseDict[basetype] if name is None else name
        return surface

    @staticmethod
    def get_dodecahedron() :
        """Numpy arrays of dodecahedron vertices (20,3) and face vertex indices (12,5)
           and edge vertex indices (30,2)
        """       
        v,f,e = SphericalSurface._get_platonic_solids_dictionary('dodeca')
        return np.array(v), np.array(f), np.array(e)

    @staticmethod
    def get_cube() :
        """Numpy arrays of cube vertices (8,3) and face vertex indices 6,5)
           and edge vertex indices (12,2)
        """
        v,f,e = SphericalSurface._get_platonic_solids_dictionary('cube')
        return np.array(v), np.array(f), np.array(e)

    @staticmethod
    def _get_platonic_solids_dictionary(getArrays = None) :
        """Constructs the dictionary of base spherical surface geometries."""
        # ....................................................
        def get_cube() :
            zeta = 1/3
            x0 = np.sqrt(2)/3
            y0 = x0*np.sqrt(3)
            chi = 2*x0
        
            cubeBaseVert = [ [0,0,1],
                             [chi,0,zeta],   [-x0,y0,zeta],  [-x0,-y0,zeta],
                             [-chi,0,-zeta], [x0,-y0,-zeta], [x0,y0,-zeta],
                             [0,0,-1] ]

            cubeBaseFaces=[ 
                [0,1,6,2],[0,2,4,3],[0,3,5,1],
                [7,6,1,5],[7,4,2,6],[7,5,3,4] ]
            
            cubeBaseEdges = [ [0,1], [0,2], [0,3], [7,4], [7,5], [7,6],
                              [1,6], [6,2], [2,4], [4,3], [3,5], [5,1] ]

            return cubeBaseVert, cubeBaseFaces, cubeBaseEdges 
        def get_dodecahedron() :
            
            zeta = 1/3
            ksi = np.sqrt(5)/3
            lam_squ = 8/9
            hEdge_squ = (1-ksi)/2

            Ax = 2/3
            Apx = 1/3
            Bx = np.sqrt(lam_squ-hEdge_squ)
            Cx = ksi 
            Dx = 1-Bx
            
            Apy = np.sqrt(3)/3
            By = np.sqrt(hEdge_squ)
            Cy = Apy
            Dy = By + Cy

            # -----------------------

            docTop = [ [0,0,1],
                [Ax, 0, ksi], [-Apx, Apy, ksi], [-Apx, -Apy, ksi],
                [Cx, Cy, zeta], [Dx, Dy, zeta], [-Bx, By, zeta],
                [-Bx, -By, zeta], [Dx, -Dy, zeta], [Cx, -Cy, zeta]
                ]
            docBottom = []
            for ver in docTop : docBottom.append( np.multiply(ver,[-1,1,-1]).tolist() )
            docbaseVert = np.concatenate((docTop,docBottom)).tolist()

            docbaseFaces = [
                [0,1,4,5,2], [0,2,6,7,3] , [0,3,8,9,1],
                [1,9,17,16,4], [2,5,15,14,6], [3,7,19,18,8],
                [10,12,16,17,13], [10,11,14,15,12], [10,13,18,19,11],
                [11,19,7,6,14],[13,17,9,8,18],[12,15,5,4,16]
            ]

            docbaseEdges = [
                [ 0, 1], [ 0, 2], [ 0, 3], [ 4, 5], [ 6, 7], [ 8, 9], 
                [ 1, 9], [ 1, 4], [ 2, 5], [ 2, 6], [ 3, 7], [ 3, 8],
                [10,11], [10,12], [10,13], [14,15], [16,17], [18,19], 
                [11,19], [11,14], [12,15], [12,16], [13,17], [13,18],
                [ 6,14], [ 4,16], [ 8,18], [ 5,15], [ 9,17], [ 7,19]  
            ]

            return docbaseVert, docbaseFaces, docbaseEdges
        def construct_tetrahedron() :
            cSqr = (np.tan(np.pi/3))*(np.tan(np.pi/3))
            zeta = (cSqr-1)/(2*cSqr)
            chi = np.sqrt(3*cSqr-1)*np.sqrt(cSqr+1)/(2*cSqr)
            x0 = chi*np.cos(np.pi/3)
            y0 = chi*np.sin(np.pi/3)
        
            baseVert = [ [0,0,1],
                         [chi,0,-zeta], [-x0,y0,-zeta], [-x0,-y0,-zeta] ]

            baseFaces=([0,1,2],[0,2,3],[0,3,1],[1,3,2])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        def construct_cube() :
            cubeBaseVert, cubeBaseFaces, cubeBaseEdges = get_cube()
            baseFaces = []
            baseVert = cubeBaseVert.copy()

            for face in cubeBaseFaces :
                vcSum = [0,0,0]
                # get square face center
                for vI in face : vcSum = np.add(vcSum,cubeBaseVert[vI])
                vcSum = np.divide(vcSum,4)
                vcCenter = vcSum/np.linalg.norm(vcSum)
                vIndex = len(baseVert)
                baseVert.append(vcCenter)
                # subdivide cube face into 4 triangles
                for i in range(0,4) :
                    newFace = []
                    newFace.append(vIndex)
                    newFace.append(face[i])
                    newFace.append(face[ (i+1)%4])
                    baseFaces.append(newFace)

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        def construct_cubeS() :
            zeta = 1/3
            x0 = np.sqrt(2)/3
            y0 = x0*np.sqrt(3)
            chi = 2*x0

            soff = 0.01
            x0S, y0S = soff/2, soff*np.sqrt(3.0)/2.0
            w = y0S/2.

            cubeBaseVert = [ [chi,-w,zeta], [chi, w,zeta],
                [-x0,y0,zeta],   [-x0,-y0,zeta],  [-chi,0,-zeta],  [x0,-y0,-zeta],   [x0,y0,-zeta],
                [ soff,   w, 1], [  x0S, y0S, 1], [ -x0S, y0S, 1], [ -soff,   0, 1], [ -x0S,-y0S, 1],
                [  x0S,-y0S, 1], [ soff,  -w, 1],
                [ soff,   w,-1], [  x0S, y0S,-1], [ -x0S, y0S,-1], [ -soff,   0,-1], [ -x0S,-y0S,-1],
                [  x0S,-y0S,-1], [ soff,  -w,-1] ]

            cubeBaseFaces=( 
                [ 7, 1, 6, 8], [ 8, 6, 2, 9], [ 9, 2, 4,10], [10, 4, 3,11], [11, 3, 5,12], [12, 5, 0,13], 
                [14,15, 6, 1], [15,16, 2, 6], [16,17, 4, 2], [17,18, 3, 4], [18,19, 5, 3], [19,20, 0, 5] )
            
            baseFaces = []
            baseVert = cubeBaseVert.copy()

            for face in cubeBaseFaces :
                vcSum = [0,0,0]
                # get square face center
                for vI in face : vcSum = np.add(vcSum,cubeBaseVert[vI])
                vcSum = np.divide(vcSum,4)
                vcCenter = vcSum/np.linalg.norm(vcSum)
                vIndex = len(baseVert)
                baseVert.append(vcCenter)
                # subdivide cube face into 4 triangles
                for i in range(0,4) :
                    newFace = []
                    newFace.append(vIndex)
                    newFace.append(face[i])
                    newFace.append(face[ (i+1)%4])
                    baseFaces.append(newFace)

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,8,1] }
        def construct_octahedron():
            
            baseVert = [ [0,0,1],
                         [1,0,0], [0,1,0], [-1,0,0], [0,-1,0],
                         [0,0,-1] ]

            baseFaces=([0,1,2],[0,2,3],[0,3,4],[0,4,1],
                [5,2,1],[5,3,2],[5,4,3],[5,1,4])

            Fo = len(baseFaces)

            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        def construct_octahedronS():
            soff = 0.01

            baseVertS = []
            x0S, y0S = soff/np.sqrt(2.0), soff/np.sqrt(2.0)

            baseVertS = [ [1,-y0S,0],
                          [1, y0S,0],       [0,1,0],        [-1,0,0],        [0,-1,0],
                          [ x0S,  y0S, 1], [-x0S,  y0S, 1], [-x0S, -y0S, 1], [ x0S, -y0S, 1],
                          [ x0S,  y0S,-1], [-x0S,  y0S,-1], [-x0S, -y0S,-1], [ x0S, -y0S,-1] ]

            baseFacesS=( [5,1,2], [6,2,3], [7,3,4], [8,4,0],
                         [2,6,5], [3,7,6], [4,8,7],
                         [9,2,1], [10,3,2], [11,4,3], [12,0,4],
                         [2,9,10], [3,10,11], [4,11,12] )

            FoS = len(baseFacesS)

            return { 'baseFaceIndices' : baseFacesS , 'baseVerticesCoor' : baseVertS, 'fevParam' : [FoS,5,1] }
        def construct_dodecahedron() :
            docbaseVert, docbaseFaces, docbaseEdges = get_dodecahedron()
            baseFaces = []
            baseVert = docbaseVert.copy()
            
            for face in docbaseFaces :
                vcSum = [0,0,0]
                # get pentagon face center
                for vI in face : vcSum = np.add(vcSum,docbaseVert[vI])
                vcSum = np.divide(vcSum,4)
                vcCenter = vcSum/np.linalg.norm(vcSum)
                vIndex = len(baseVert)
                baseVert.append(vcCenter)
                # subdivide pentagon face into 5 triangles
                for i in range(0,5) :
                    newFace = []
                    newFace.append(vIndex)
                    newFace.append(face[i])
                    newFace.append(face[ (i+1)%5])
                    baseFaces.append(newFace)

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        def construct_icosahedron():
            cSqr = (np.tan(np.pi/5))*(np.tan(np.pi/5))
            zeta = (1-cSqr)/(2*cSqr)
            chi = np.sqrt(3*cSqr-1)*np.sqrt(cSqr+1)/(2*cSqr)
            thetaInc = 2*np.pi/5
            thetaOff = thetaInc/2

            baseVert = []
            baseVert.append([0,0,1])
            for i in range(0,5) :
                theta = i*thetaInc
                baseVert.append( [ chi*np.cos(theta), chi*np.sin(theta), zeta] )
            for i in range(0,5) :
                theta = thetaOff+i*thetaInc
                baseVert.append( [ chi*np.cos(theta), chi*np.sin(theta), -zeta] )
            baseVert.append([0,0,-1])

            baseFaces=([0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,1],
                [6,2,1],[7,3,2],[8,4,3],[9,5,4],[10,1,5],
                [2,6,7],[3,7,8],[4,8,9],[5,9,10],[1,10,6],
                [11,7,6],[11,8,7],[11,9,8],[11,10,9],[11,6,10])

            Fo = len(baseFaces)
            return { 'baseFaceIndices' : baseFaces , 'baseVerticesCoor' : baseVert, 'fevParam' : [Fo,0,2] }
        # ....................................................
        if getArrays == 'cube' : return get_cube()
        if getArrays == 'dodeca' : return get_dodecahedron()

        surfacesDictionary = {}
        surfacesDictionary['tetra']  = construct_tetrahedron()
        surfacesDictionary['cube']   = construct_cube()
        surfacesDictionary['cube_s'] = construct_cubeS()
        surfacesDictionary['octa']   = construct_octahedron()
        surfacesDictionary['octa_s'] = construct_octahedronS()
        surfacesDictionary['dodeca'] = construct_dodecahedron()
        surfacesDictionary['icosa']  = construct_icosahedron()
        # ... 'special case' only for fev method for spherical surfaces.
        surfacesDictionary['octa_c'] = { 'fevParam' : [16,5,1] }
        surfacesDictionary['cube_c'] = { 'fevParam' : [12,4,1] }
        # ... 'special case' only for fev method for platonic surfaces.
        surfacesDictionary['cube_a'] = { 'fevParam' : [12,0,2] }
        surfacesDictionary['cube*'] = { 'fevParam' : [12,0,2] }
        surfacesDictionary['dodeca*'] = { 'fevParam' : [36,0,2] }

        return surfacesDictionary

    @staticmethod
    def _midVectorFun(vectA, vectB) :
        """Mid-point between two points on a unit sphere."""
        mid = np.add(vectA,vectB)
        return mid/np.linalg.norm(mid)
    
    @staticmethod
    def _map3Dto2Dsquare(xyz) :
        """Override map surface to rectagular unit coordinates (0,1)."""
        x,y,z = xyz
        t = np.arctan2(y,x)
        theta = np.where( t<0, 2.0*np.pi + t ,t)
        phi = np.arcsin(z)
        a = theta/(2*np.pi)
        b =  phi/np.pi + 0.5
        a = np.clip(a,0.0,1.0) 
        b = np.clip(b,0.0,1.0) 
        return [a,b]

    @staticmethod
    def coor_convert(xyz, tocart=False) :
        """
        Transformation between spherical and Cartesian coordinates.
        The azimuthal coordinate is in radians with domain 0 to 2pi.
        The polar angle is in radians with domain 0 to pi.

        Parameters
        ----------
        xyz : 3xN array
            N number of vectors in either spherical or cartesian coordinates,

        tocart : { True, False }, default: False
            If True, input is spherical and output is cartesian.
            If False, input is cartesian and output is spherical.

        Returns
        -------
        list: 3xN array
            Transformed coordinates. 

        """

        if tocart :
            r, theta, phi = xyz
            x = r*np.sin(phi)*np.cos(theta)
            y = r*np.sin(phi)*np.sin(theta)
            z = r*np.cos(phi)
            abc = [x,y,z]
        else:
            x,y,z = xyz
            R = np.sqrt( x*x + y*y + z*z )
            phi = np.nan_to_num(np.arccos(z/R))
            theta = np.arctan2(y,x)
            theta = np.where(theta<0,theta+2*np.pi,theta)
            abc = [R,theta,phi]
        abc = np.array(abc)
        return abc

    _base_surfaces = _get_platonic_solids_dictionary.__func__()
    _default_base = 'icosa'  # default assignment is indicated in __init Docstring.
    _geomType = _COORSYS['SPHERICAL']

    @classmethod
    def _uvw_to_xyz(cls, uvw, rtp) :
        """Override rotational transformation at a spherical coordinate.""" 
        return super()._disp_to_xyz( uvw, rtp[1], rtp[2] )

    @classmethod
    def grid(cls,nlat, nlng, basetype='d', minrad=None, name=None, **kwargs) :
        """
        SphericalSurface with latitude and longitude faces.

        Parameters
        ----------
        nlat : integer, required.
            Number of latitude subdivisions.  Minimum value is 2

        nlng : integer, required.
            Number of longitude subdivisions. Minimum value is 3.

        basetype : {'d','s','q','w','r','x'}, optional, default: 'd'
            Starting surface geometries from which the surface is
            constructed indicating face shape and if split.
        
        minrad : scalar, optional, default: 0.01
            The minimum distance from the z-axis of any vertex coordinate.
            Use is dependent on basetype.

        name : string, optional, default: None.
            Descriptive identifier.

        Raises
        ------
        ValueError
            If nlat is not an integer greater or equal to 2.
            If nlng is not an integer greater or equal to 3.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.
        
        Returns
        -------
        surface : SphericalSurface object.

        """
        Nfaces = { 'd': 6, 's': 6, 'q':12, 'w':12, 'r':6, 'x':6  }  # 2,3
        splitSize = None
        nrad, nang = nlat, nlng
        surface = cls._square_net(nrad, nang, basetype,minrad, splitSize, name, **kwargs)
        surface._postProc_surfaceColors(True)
        surface._rez = surface._calc_rez( 20 )
        return surface

    @classmethod
    def pntsurf(cls, rtp_points, name=None, **kargs) :
        """
        SphericalSurface from points in spherical coordinates.

        Parameters
        ----------
        rtp_points : array, required.
            An N X 3, or N X 4, array of N spherical coordinate
            points. For an N X 4 array, the 4th is the vertex
            scalar value.

        name : string, optional, default: None.
            Descriptive identifier for the point surface.

        Returns
        -------
        surface : SphericalSurface object

        """
        dataPts = np.array(rtp_points)
        # update for v_1.3.0
        dataPts, vVal = Surface3DCollection._extract_vals_from_inputCoor(dataPts)

        sphPts = dataPts.T
        xyzPts = SphericalSurface.coor_convert(dataPts.T,True)
        rmax = np.amax(sphPts[0])
        sphPts[0] = rmax
        xyz_sphPts = SphericalSurface.coor_convert(sphPts,True)
        hull_surface = Surface3DCollection.chull(xyz_sphPts.T)
        verts = xyzPts.T
        faceIndices = hull_surface.fvIndices

        if name is None : name = 'point surface'
        surface = SphericalSurface(name=name,**kargs) 
        surface.vertexCoor = verts
        surface.fvIndices = faceIndices
        surface.evIndices = surface._get_edges_from_faces(faceIndices)
        surface.vfIndicesList = surface._vertexCommonFaceIndices()
        surface.efIndicesList = surface._edgeCommonFaceIndices()
        surface.set_verts( verts[faceIndices] )
        surface._set_geometric_bounds()
        surface._rez = surface._calc_rez(20)
        surface._basetype = 'pntsurf'
        if vVal is not None : surface.set_vertvals(vVal)   # update for v_1.3.0
        return surface

    @classmethod
    def _Xfev(cls, rez, basetype=None, minrad=None) :
        """
        Calculates number of faces, edges, and vertices of a SphericalSurface.

        Parameters
        ----------
        rez : integer
            Number of recursive subdivisions of a triangulated base faces.
        
        basetype : string
            Starting surface geometries from which the surface would be
            constructed using recursive subdivisions of the triangular faces.
        
        Returns
        -------
        list: the number of faces, edges, and vertices.

        """
        if isinstance(rez,(list,tuple)) : # check for globe...
            nlat = rez[0]
            nlng = rez[1]
            # only check is ending in s
            isSplit = basetype[ len(basetype)-1] == 's'
            if minrad is None :
                f = 2*nlat*nlng - 2*nlng
                e = 3*nlat*nlng - 3*nlng
                v = nlat*nlng - nlng + 2
                if isSplit :
                    e += nlat
                    v += nlat - 1
            else :  # any value will trigger this condition.
                f = 2*nlat*nlng
                e = 3*nlat*nlng + nlng
                v = nlat*nlng + nlng
                if isSplit :
                    e += nlat    
                    v += nlat + 1    
            return [f,e,v]

        return super()._fev(rez,basetype,cls)

    @staticmethod
    def _spherical_split_cylinder(rez,basetype,minrad, **kwargs) :
        """
        Transform cylindrical base surface to spherical coordinates. 

        note :  This method provides an alternative to the 
                spherical basetypes of 'hex_s' and 'square_s'.
                The grids constructed here may provide smoother
                geometries depending on the mapping functions.
        """

        # .................................................
        def cyl2sph(rtz):
            r,t,z = rtz
            minphi = np.arccos(minrad)
            phi = z*minphi
            Z = np.sin(phi)
            R = np.cos(phi)
            return R,t,Z
        # .................................................
        if basetype == 'octa_c' : cylbase = 'squ_s' 
        if basetype == 'cube_c' :    cylbase = 'tri_s'
        cyl_obj = CylindricalSurface(rez,cylbase, **kwargs).map_geom_from_op(cyl2sph)
        baseFaceVertexIndices = cyl_obj._base_surfaces[cylbase] ['baseFaceIndices']    
        fvIndices = cyl_obj.fvIndices
        vertexCoor = cyl_obj.vertexCoor
        evIndices = cyl_obj.evIndices
        return vertexCoor, fvIndices, evIndices, baseFaceVertexIndices,

    def __init__(self, rez=0, basetype=None, minrad=None, name=None, **kwargs):
        """
        Create spherical surface of unit radius.

        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of the triangulated base faces.
            Rez values range from 0 to 7.
            
        basetype : {'tetra','octa','icosa','cube','dodeca', \
            'octa_s','cube_s','octa_c','cube_c'}, optional, default: 'icosa'
            Starting surface geometries from which the surface is
            constructed using recursive subdivisions of the triangular faces.
            All basetypes are based on the five platonic solids.

        minrad : scalar, optional, default: 0.01
            For basetypes 'octa_c' and 'cube_c, the minimum
            distance from the z-axis of any vertex coordinate.

        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Raises
        ------
        ValueError
            If rez is not an integer in range 0 to 7.
            If basetype is not recognized for this surface constructor.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.Surface3DCollection'.

        """

        if minrad is None : minrad = _MINRAD
        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))
        if basetype is None : basetype = self._default_base
        # ... check for 'special case' (ugly code but needed, MROI)
        if basetype == 'octa_c' or basetype == 'cube_c' :
            vertexCoor, faceIndices, edgeIndices, bfvi = self._spherical_split_cylinder(rez,basetype,minrad, **kwargs)
            baseFaceVertexIndices = bfvi
        else :
            if basetype not in self._base_surfaces :
                raise ValueError('Basetype {} is not recognized. Possible values are: {}'.format(basetype,list(self._base_surfaces)))
            baseSurfObj = self._base_surfaces[basetype]
            baseVcoor = baseSurfObj['baseVerticesCoor']
            baseFaceVertexIndices =baseSurfObj['baseFaceIndices']
            indexObj,vertexCoor = self._init_triangulateBase(rez,baseVcoor,baseFaceVertexIndices, self._midVectorFun)
            faceIndices = indexObj['face']
            edgeIndices = indexObj['edge']
        
        super().__init__(vertexCoor, faceIndices, name, **kwargs)

        self.coorType = _COORSYS["SPHERICAL"]
        self._normal_scale = _normLength( len(faceIndices), 0.61, 1.52)
        self._rez = rez
        self._basetype = basetype
        self.disp_Vector = np.array( [1,1,1] )
        return None
    
    def domain(self, radius=None) :
        """
        Set the domain of the spherical surface.

        Used for setting the radius of the base sphere.

        Parameters
        ----------
        radius: number, default : 1
            Radius of the spherical surface.
            If set to None, the radius is unchanged.

        Returns
        -------
        self : SphericalSurface object

        """
        def setDomain(rtp, mr ) :
            r,t,p = rtp
            R = mr*r
            return R,t,p       
        def getCont(lims) :
            m = 1.0
            if lims is not None : m = lims
            return m
        mr = getCont(radius)
        self.map_geom_from_op(lambda A : setDomain(A, mr ))
        return self



class Vector3DCollection(Line3DCollection) :
    """
    Collection of 3D vectors represented as arrows.

    """

    # FutDev: additional methods to add functionality.

    def __init__(self, location, vect , alr=None, name=None, **kwargs) :
        """
        Create a collection of 3D vectors.

        Parameters
        ----------
        location : N x 3 float array
            Cartesian coordinate location (tails) for N number of vectors.

        vect : N x 3 float array
            N number of vectors in Cartesian coordinates.

        alr : scalar, optional, default: 0.25
            Axis length ratio, head size to vector magnitude.

        name : string identifier

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to mpl_toolkits.mplot3d.art3d.Line3DCollection.
            Valid keywords include: colors, linewidths.
        
        """
        # note: self coor types primarially used for export info.

        self.coorType = _COORSYS["XYZ"]
        self.uvwOrientationType = _COORSYS["XYZ"]

        if alr is None : alr = _DFT_ALR
        self._alr = alr

        self.baseLocation = np.array(location, dtype=float)
        self.baseVect = np.array(vect, dtype=float)

        lines,XYZ,UVW = self._quiverLines(location,vect,alr)
        # note: _quiverLines may remove zero length vect
        self.location = XYZ
        self.vect = UVW

        if 'color' in kwargs :
            try:
                test = colors.to_rgba( kwargs['color'] )
            except:
                del kwargs['color']
                warnings.warn("Only a SINGLE color value may be used for vector instantiation")

        super().__init__(lines, **kwargs)
        #self._geomName = name    # Redunant from next line, internal tracking to determine if set...
        self.name = name
        self.valuesName = None
        self.baseColor = np.array(self.get_color()[0]).flatten().tolist()


        self._bounds = {}
        _set_geomBounds(self.location,self._bounds)
        self._bounds['vlim'] =     _DFT_VLIM         #.. from operation
        self._bounds['vertvlim'] = _DFT_VERTVLIM     #.. from direction or magnitude

    def _set_vectColor(self, color) :
        
        # note: for a list of four colors [ A,B,C,D ], result is
        #       [ A,B,C,D, A,A, B,B, C,C, D,D ].
        #       line segment colors, then a pair for each arrow head (2 lines).

        head = np.arange(len(color))
        tail = np.ravel( [ [head],[head] ] , order='F')
        index = np.append(head,tail)
        #self._edgecolors = color[index]
        self.set_color(color[index])     # update for v_1.2.0
        return

    def __str__(self) :
        # note: each vector representation is composed of 3 lines
        numbVect = int(len(self._segments3d)/3)
        name = self.__class__.__name__
        val = name + ':  vectors: {}'.format(numbVect)
        return val

    def _quiverLines( self, vC, vN, alr) :
        # ..................................................
        def quiver( *args, length=1, arrow_length_ratio=alr, pivot='tail', normalize=False):
            # taken directly from mpl_toolkits.mplot3d.axes3d source 
            # and removed axis reference.  Acknowledgement:
            # Copyright (c) 2012-2013 Matplotlib Development Team; All Rights Reserved
            # ...............................................
            def calc_arrow(uvw, angle=15):
                """
                To calculate the arrow head. uvw should be a unit vector.
                We normalize it here:
                """
                # get unit direction vector perpendicular to (u,v,w)
                norm = np.linalg.norm(uvw[:2])
                if norm > 0:
                    x = uvw[1] / norm
                    y = -uvw[0] / norm
                else:
                    x, y = 0, 1

                # compute the two arrowhead direction unit vectors
                ra = np.radians(angle)
                c = np.cos(ra)
                s = np.sin(ra)

                # construct the rotation matrices
                Rpos = np.array([[c+(x**2)*(1-c), x*y*(1-c), y*s],
                                [y*x*(1-c), c+(y**2)*(1-c), -x*s],
                                [-y*s, x*s, c]])
                # opposite rotation negates all the sin terms
                Rneg = Rpos.copy()
                Rneg[[0, 1, 2, 2], [2, 2, 0, 1]] = \
                    -Rneg[[0, 1, 2, 2], [2, 2, 0, 1]]

                # multiply them to get the rotated vector
                return Rpos.dot(uvw), Rneg.dot(uvw)
            # ...............................................
            argi = 6
            if len(args) < argi:
                raise ValueError('Wrong number of arguments. Expected %d got %d' %
                                (argi, len(args)))
            input_args = args[:argi]
            masks = [k.mask for k in input_args
                    if isinstance(k, np.ma.MaskedArray)]
            bcast = np.broadcast_arrays(*input_args, *masks)
            input_args = bcast[:argi]
            masks = bcast[argi:]
            if masks:
                mask = reduce(np.logical_or, masks)
                input_args = [np.ma.array(k, mask=mask).compressed()
                            for k in input_args]
            else:
                input_args = [np.ravel(k) for k in input_args]

            if any(len(v) == 0 for v in input_args):
                return []
            shaft_dt = np.array([0, length])
            arrow_dt = shaft_dt * arrow_length_ratio
            if pivot == 'tail':
                shaft_dt -= length
            elif pivot == 'middle':
                shaft_dt -= length/2.
            elif pivot != 'tip':
                raise ValueError('Invalid pivot argument: ' + str(pivot))
            XYZ = np.column_stack(input_args[:3])
            UVW = np.column_stack(input_args[3:argi]).astype(float)
            norm = np.linalg.norm(UVW, axis=1)
            mask = norm > 0
            XYZ = XYZ[mask]
            if normalize:
                UVW = UVW[mask] / norm[mask].reshape((-1, 1))
            else:
                UVW = UVW[mask]
            if len(XYZ) > 0:
                shafts = (XYZ - np.multiply.outer(shaft_dt, UVW)).swapaxes(0, 1)
                head_dirs = np.array([calc_arrow(d) for d in UVW])
                heads = shafts[:, :1] - np.multiply.outer(arrow_dt, head_dirs)
                heads.shape = (len(arrow_dt), -1, 3)
                heads = heads.swapaxes(0, 1)
                lines = [*shafts, *heads]
            else:
                lines = []
            return lines , XYZ,UVW   # <--- added return of XYZ,UVW
        # ..................................................
        x,y,z = np.transpose( vC )
        u,v,w = np.transpose( vN )
        lines, XYZ, UVW = quiver(x,y,z,u,v,w)
        return lines, XYZ, UVW

    @property
    def alr(self) :
        '''Arrow head to length ratio.'''
        return self._alr

    @alr.setter
    def alr(self, val) :
        location = self.baseLocation
        vect = self.baseVect
        if val is None : val = _DFT_ALR
        self._alr = val
        lines,_,_ = self._quiverLines(location,vect,self._alr)
        self.set_segments(lines)
        return

    @property
    def vlim(self) :
        '''Range of values associated with color'''
        return self._bounds['vlim']

    @vlim.setter
    def vlim(self,val) :
        self._bounds['vlim'] = val
        return

    @property
    def name(self) :
        '''Descriptive identifier for the vector geometry.'''
        if self._geomName is None : return ''
        return self._geomName

    @name.setter
    def name(self, val) :
        self._geomName = val
        self.set_label( '' if val is None else val)
        return

    @property
    def cname(self) :
        '''Descriptive identifier for values indicated by color.'''
        if self.valuesName is None : return ''
        return self.valuesName

    @cname.setter
    def cname(self, val) :
        self.valuesName = val
        return

    @property
    def bounds(self) :
        """
        Dictionary of vector geometric and value ranges.
        
        Each dictionary value is a 2 float array of minimum
        and maximum values of the vector location.  Keys are:

        'xlim' : x-coordinate

        'ylim' : y-coordinate

        'zlim' : z-coordinate

        'r_xy' : radial distance from the z axis

        'rorg' : radial distance from the origin

        'vlim' : value functional assignments.

        'vertvlim: magnitude
            
        Values are assigned from the geometry and color mapping methods.

        """

        return self._bounds

    def map_color_from_op( self, operation, rgb=True, cname=None ) :
        """
        Assignment of vector color from a function.

        Vector colors are assigned from a function
        of direction and location coordinates.

        Parameters
        ----------
        operation : function object
            Function that takes two arguments, both
            a 3xN Numpy array of xyz coordinates.
            The first and second arguments are the location
            and direction, respectively.
            The function returns a 3xN color value.

        rgb : bool {True, False}, optional, default: True
            By default, RGB color values are returned by the
            operation function.  If set False, the operation
            returns HSV color values.

        Returns
        -------
        self : Vector3DCollection object

        """
        xyz = np.transpose(self.location)
        uvw = np.transpose(self.vect)
        # .....
        colors = np.array(operation(xyz,uvw))
        colors = np.transpose(colors)
        colors = np.clip(colors,0,1)
        if rgb : RGB = colors
        else :   RGB = cm.colors.hsv_to_rgb(colors)

        if RGB.shape[1] == 3 :
            ones = np.ones(RGB.shape[0])[:,np.newaxis]
            RGB = np.concatenate((RGB,ones),axis=1)

        self._set_vectColor(RGB)
        cname = _getFunctionName(operation,cname)
        if cname is not None : self.valuesName = cname
     
        return self

    def map_cmap_from_op( self, operation, cmap=None, cname=None ) :
        """
        Functional assignment of a vector color from a color map.

        Location and direction coordinates are used to calculate a scalar 
        which is then used to assign face colors from a
        colormap.

        Parameters
        ----------
        operation : function object
            Function that takes two arguments, both
            a 3xN Numpy array of xyz coordinates.
            The first and second arguments are the location
            and direction, respectively.
            The function returns a Numpy array of scalar values.

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the function return values to colors.
        
        Returns
        -------
        self : Vector3DCollection object

        """
        
        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cm_colorMap = self.get_cmap()

        xyz = np.transpose(self.location)
        uvw = np.transpose(self.vect)
        # .....
        v = np.array(operation(xyz,uvw))
        norm = colors.Normalize(vmin=v.min(),vmax=v.max())
        color = cm_colorMap(norm(v))
        self._set_vectColor(color)
        cname = _getFunctionName(operation,cname)
        if cname is not None: self.valuesName = cname
        self._bounds['vlim'] = [ v.min(), v.max() ]

        return self

    def map_cmap_from_magnitude( self, cmap=None, cname=None ) :
        """
        Vector color assignment using vector magnitude.

        Parameters
        ----------
        cmap : str or Colormap, optional, default: 'viridis'
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the dot product values to colors.
        
        Returns
        -------
        self : Vector3DCollection object

        """

        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cm_colorMap = self.get_cmap()

        mag = np.linalg.norm(self.vect,axis=1)
        norm = colors.Normalize(vmin=mag.min(),vmax=mag.max())
        color = cm_colorMap(norm(mag))
        self._set_vectColor(color)
        self.cname = cname
        if cname is None: self.cname = 'Magnitude'
        self._bounds['vertvlim'] = [ mag.min(), mag.max() ]

        return self
    
    def map_cmap_from_direction( self, cmap=None, direction=[1,1,1], cname=None ) :
        """
        Vector color assignment using vector direction relative to direction argument.

        The dot product of vector direction with the argument direction 
        is used to assign vector colors from a colormap.

        Parameters
        ----------
        cmap : str or Colormap, optional, default: 'viridis'
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the dot product values to colors.
        
        direction : list of size 3, optional, default: [1,1,1]
            A 3D vector in xyz Cartesian coordinates designating
            the reference direction.

        refCoor : string, optional, default: "XYZ"
            Direction coordinate system for the evaluation.
            (Not implimented)

        Returns
        -------
        self : Vector3DCollection object

        """

        # FutDev:  allow input direction to be in other coor systems (refCoor)
        #.....................................    
        def getColorMap(fvo,cmap, direction) :
            unitVector = lambda v : np.divide( v, np.linalg.norm(v) )
            incidentLight = unitVector( direction )
            d = np.dot(fvo,incidentLight)
            v = np.add(d,1)
            v = np.divide(v,2)
            v = np.abs(v)
            norm = colors.Normalize(vmin=v.min(),vmax=v.max())
            colorMap = cmap(norm(v))
            self._bounds['vertvlim'] = [ v.min(), v.max() ]
            return colorMap
        #.....................................

        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cm_colorMap = self.get_cmap()

        color = getColorMap(self.vect, cm_colorMap, direction)
        self._set_vectColor(color)
        self.cname = cname
        if cname is None: self.cname = 'Direction, ' + str(direction)

        return self



class ColorLine3DCollection(Line3DCollection) :
    """
    Base class for non-uniformly colored 3D lines.

    """

    # FutDev: This class leaves a lot of room for improvement.
    # Need better use of Numpy and better approach to
    # differentiating between collections, lines & segments
    # to accomodate Matplot rendering and future export. MROI

    @staticmethod
    def meshgrid(X,Y,Z, name=None) :
        """
        ColorLine3DCollection, a set of xy slice lines based on a meshgrid.
        
        Parameters
        ----------
        X,Y,Z : N x M arrays
            Cartesian x,y,z coordinate values.

        name : string, optional, default: None.
            Descriptive identifier for the mesh surface.

        Returns
        -------
        surface : ColorLine3DCollection object
        
        """

        # Needs cleanup MROI
        x_lineVertsPts = np.array(X)[0]
        y_lineVertsPts = np.array(Y)[:,0]
        z_lines = np.array(Z)
        numb_X_lines = len(x_lineVertsPts)
        numb_Y_lines = len(y_lineVertsPts)
        indices = [None]*( numb_X_lines + numb_Y_lines )
        index = 0

        verts = []
        # ..................................................
        numVerts = numb_Y_lines
        lineVertsPts = y_lineVertsPts
        z_X_lines = z_lines
        for iX in range(numb_X_lines) :
            for iY in range(numVerts) :
                verts.append( [ x_lineVertsPts[iX] , lineVertsPts[iY], z_X_lines[iY,iX] ])

        for i in range(0,numb_X_lines) :
            indices[i] = [None]*(numVerts)
            for j in range(numVerts) :
                indices[i][j] = index
                index += 1
        # ..................................................
        numVerts = numb_X_lines
        lineVertsPts = x_lineVertsPts
        z_Y_lines = z_lines
        for iY in range(numb_Y_lines) :
            for iX in range(numVerts) :
                verts.append( [  lineVertsPts[iX], y_lineVertsPts[iY] ,z_Y_lines[iY,iX] ])

        for i in range(numb_X_lines,numb_X_lines+numb_Y_lines) :
            indices[i] = [None]*(numVerts)
            for j in range(numVerts) :
                indices[i][j] = index
                index += 1

        verts = np.array(verts)

        return ColorLine3DCollection(verts,indices,name=name)

    def __init__(self, vertexCoor, segmIndices, name=None, **kwargs) :
        """
        A set 3D lines, each line with a varying number of sequential vertices.

        Parameters
        ----------
        vertexCoor : V x 3 float array-like 
            An array of V number of xyz vertex coordinates.

        segmIndices : L x Nl int list
            An list of L number of Nl line vertex indices, where
            Nl is the number of vertices for the Lth line.

        name : string, optional, default: None.
            Descriptive identifier.

        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to mpl_toolkits.mplot3d.art3d.Line3DCollection.
            Valid keywords include: colors, linewidths.
        
        """

        self.vertexCoor = np.array(vertexCoor,dtype=float)

        if segmIndices is None :
            self.lineIndices = [ [*range(len(vertexCoor))] ]

            linesegs = self.lineIndices
            indices = []
            for line in linesegs :
                numSegs = len(line) - 1
                for seg in range (numSegs) :
                    t = [ line[seg],line[seg+1] ]
                    indices.append(t)
            self.segmIndices = indices

        else : 
            # all indices must be lists, not numpy arrays.
            segmIndices = list(segmIndices)

            # segmIndices will be shredded for any color operation
            # forming multiple single segment lines.
            self.segmIndices = segmIndices.copy()

            # lineIndices are appended for an add operation,
            # where lineIndices preserve individual 'lines' even though
            # all lines are broken into single segment lines omce shredded.
            self.lineIndices = segmIndices.copy()

        self.isShredded = False

        S = self.segmIndices
        v = self.vertexCoor

        segments = [None]*len(S)
        for verlist in range(len(S)) :
            verArr = [None]*len( S[verlist] )
            for verInx in range( len( S[verlist] ) ) :
                verArr[verInx] = v[ S[verlist][verInx] ]
            segments[verlist] = verArr

        super().__init__(segments, **kwargs)
        #self._geomName = None    # Redundant from next line, internal tracking to determine if set...
        self.name = name
        self.valuesName = None

        # ...............................................................
        self.coorType = _COORSYS["XYZ"]
        self.baseColor = np.array(self.get_color()[0]).flatten().tolist()
        self.vertexColor = np.array( [ self.baseColor ] )
        # following assigned when surface color applied using an operation
        self.vertexValues = None
        self.segValues = None
        # ...............................................................
        self.set_joinstyle('round')
        self.set_capstyle('round')
        self.scale = [1,1,1]                      # for axis when transformed.
        self.translation = [0,0,0]                # for axis when transformed.
        self.rotation =  np.identity(3).tolist()  # for axis when transformed.
        self.restoredClip = None

        self.rez = None
        self._bounds = {}
        self._set_geometric_bounds()
        self._bounds['vlim'] =     _DFT_VLIM 
        self._bounds['vertvlim'] = _DFT_VERTVLIM
        self._planeNormal = None  #.. set for planar contours, for FutDev

    @staticmethod
    def _midVectorFun(vectA, vectB) :
        """Mid-point between two points in the xy plane."""
        mid = np.add(vectA,vectB)
        mid = np.multiply(0.5,mid)
        return mid

    @staticmethod
    def _extractLines(lineSegments) :
        # Edge array sorting to line arrays, which may be loops. Edge arrays
        # must have two indices in consistent order.
        # This is an 'unshred' function.
        # ..................................................................
        def extractLine(edgeDict) :
            # use a dictionary to build a 'sort-of' single linked list
            # note: index-list would conain many emply values.
            edgeDict_fl = edgeDict
            k = list(edgeDict.keys())[0]
            v = list(edgeDict.values())[0]
            line = [k,v]
            edgeDict_fl.pop(k)

            # add segments to the tail :
            continueBuild = True
            while continueBuild :
                continueBuild = False
                end = line[len(line)-1]
                if end in edgeDict_fl :
                    line.append(edgeDict_fl[end])
                    edgeDict_fl.pop(end)
                    continueBuild = True
            if line[len(line)-1] == line[0] : 
                return line,edgeDict_fl
            # then line is not a loop so continue, adding segments to head :
            
            edgeDict_lf = { value:key for key,value in edgeDict_fl.items()}
            continueBuild = True
            while continueBuild :
                continueBuild = False
                start = line[0]
                if start in edgeDict_lf :
                    line.insert(0,edgeDict_lf[start])
                    edgeDict_lf.pop(start)
                    continueBuild = True
            edgeDict_nxt = { value:key for key,value in edgeDict_lf.items() }        
            return line, edgeDict_nxt
        # ..................................................................
        edgeDict =  { litem[0]:litem[1]  for litem in lineSegments }
        lines = []
        while(len(edgeDict)>0) :
            line,edgeDict = extractLine(edgeDict)
            lines.append(line)
        return lines

    def __add__(self, othr) :
        """
        Combine two line objects into a single line object.

        Parameters
        ----------
        othr : ColorLine3DCollection object 
        
        Returns
        -------
        ColorLine3DCollection object

        """

        # ............................................................................
        def get_lineColors( colors, lines) :
            colors = np.array(colors)
            nlines = len(lines)
            if len(colors) != nlines :
                colors = np.tile(colors,(nlines,1))
            return colors
        # ............................................................................
        
        if not isinstance(othr, ColorLine3DCollection) :
            raise ValueError('Add operations can only apply between ColorLine3DCollection objects.')

        # vertexCoor are array-like.
        topVerCoor = self.vertexCoor
        botVerCoor = othr.vertexCoor
        totVerCoor = np.append(topVerCoor,botVerCoor,axis=0)
        nextIndex = len(topVerCoor)

        # lineIndices are int Lists
        toplineIndices = self.lineIndices
        botlineIndices = copy.copy(othr.lineIndices)
        for i in range(len(botlineIndices)) : 
            botlineIndices[i] = np.add( botlineIndices[i] , nextIndex ).tolist()
        totLineIndices = copy.copy(toplineIndices)
        for line in botlineIndices :
            totLineIndices.append(line)
        
        # segmIndices are int Lists
        topSegmIndices = self.segmIndices
        botSegmIndices = copy.copy(othr.segmIndices)
        for i in range(len(botSegmIndices)) : 
            botSegmIndices[i] = np.add( botSegmIndices[i] , nextIndex ).tolist()
        totSegmIndices = copy.copy(topSegmIndices)
        for line in botSegmIndices :
            totSegmIndices.append(line)

        # edgecolors are array-like
        topEdgeColors = get_lineColors(self._ledgecolors,topSegmIndices)  # update for v_1.2.0  
        botEdgeColors = get_lineColors(othr._ledgecolors,botSegmIndices)  # update for v_1.2.0
        totEdgeColors = np.append(topEdgeColors,botEdgeColors,axis=0)

        # FutDev: 'add' will not retain line collections.  Needed for export
        #          but sufficient for rendering. MROI
        newLine = ColorLine3DCollection(totVerCoor,totSegmIndices,color=totEdgeColors)
    
        newLine.isShredded = self.isShredded and othr.isShredded
        newLine._planeNormal = None  #.. could check here if same normals.
        newLine._set_geometric_bounds()

        return newLine

    def __str__(self) :
        # note: each line representation is composed of multiple segment lines.
        numVerts = len(self.vertexCoor)
        numLines = len(self.lineIndices)
        numSeg = 0
        for line in self._segments3d :
            numSeg += len(line) -1
        name = self.__class__.__name__
        if self.rez is not None : name = name + ' ({}):'.format(self.rez)
        val = name + ': lines: {}, segs: {}, verts: {}'.format(numLines,numSeg,numVerts)
        return val

    def _postProc_segmentColors(self) :
        """ Matplotlib will 'fillin' colors during rendering using the
            facecolor array. (this is NOT the assigned color to edges).
            However, S3Dlib will use the inherited _edgecolors 
            for the assigned edeg colors.
            Whenever color is 'externally' assigned or reassigned,
            need to directly assign colors for each edge for any
            method that uses edgecolors.  This occurs initially
            in the method before colors are assigned            
            (set_color, +, clipping, shading, edges, etc.)
            
        """
        N = len(self.segmIndices)    # number of segments
        initEC = self._ledgecolors   # update for v_1.2.0
        inLen = len(initEC)          # number of edgecolors
        if inLen == N : return
        if isinstance(initEC,np.ndarray) : initEC = initEC.tolist()
        n = math.ceil(N/inLen)
        colorcoll = initEC*n
      
        self.set_color(colorcoll[0:N])
        return

    def _shred(self) :
        # convert multi-segment lines to multiple one-segment lines.
        # Required for application of non-uniform line color.
        
        linesegs = self.lineIndices
        indices = []
        for line in linesegs :
            numSegs = len(line) - 1
            for seg in range (numSegs) :
                t = [ line[seg],line[seg+1] ]
                indices.append(t)
        segments = self.vertexCoor[np.array(indices)]
        self.segmIndices = indices
        self.set_segments(segments)
        self.isShredded = True
        return segments

    def _shred_segm_color(self,numDiv) :
        initColors = self._ledgecolors  # colors before line shredding  # update for v_1.2.0
        numColors = len(initColors)
        linesegs = self.segmIndices
        colors = []
        colorIndex = 0
        for line in linesegs :
            currentColor = initColors[colorIndex]
            for seg in range (numDiv) :
                colors.append(currentColor)
            colorIndex = (colorIndex+1)%numColors
        colors = np.array(colors)     # colors after line shredding

        return colors

    def _shred_color(self) :
        # Since line colors are in a SINGLE List of colors
        # for mulitple lines with multiple segment colors,
        # colors for a multi-segment line of a single color
        # is reconstructed into multiple colors of the
        # same color.  ( see _shred )
        # Individual segemnt colors are needed for line
        # shading, clipping, and alpha setting.

        initColors = self._ledgecolors  # colors before line shredding  # update for v_1.2.0
        numColors = len(initColors)
        linesegs = self.lineIndices
        colors = []
        colorIndex = 0
        for line in linesegs :
            numSegs = len(line) - 1
            currentColor = initColors[colorIndex]
            for seg in range (numSegs) :
                colors.append(currentColor)
            colorIndex = (colorIndex+1)%numColors
        colors = np.array(colors)     # colors after line shredding

        return colors

    def _get_segment_centers(self,segments) :
        # segment centers from segment point coordinates.
        # Assumed segments of a shredded line.
        vfT = np.transpose(segments, (1,2,0) )
        numVerts = vfT.shape[0]
        sumV = np.sum(vfT, axis=0)
        return np.transpose(sumV)/numVerts

    def _get_segment_directions(self,segments) :
        # segment unit directions from segment point coordinates.
        # Assumed segments of a shredded line.
        vfT = np.transpose(segments, (1,0,2) )
        direction = vfT[1] - vfT[0]
        sz = np.linalg.norm(direction, axis=1)[:,np.newaxis]
        unitVector = np.divide( direction, sz )
        return unitVector

    def _set_geometric_bounds(self) :
        _set_geomBounds(self.vertexCoor,self._bounds)
        return    
    
    def _get_segments( self, v, S  ) :
        # Since line segments are in a List of multiple lines,
        # line vertices are set per line in the line collection
        # to reconstruct the segment coordinates ( self._segments3d )

        # DevNote: got to be a better way of doing this, but line
        #          sets are lists with varying number of segments
        #          until shredded, then all the same number of segs.

        if self.isShredded :  return v[np.array(S)]

        segments = [None]*len(S)
        for verlist in range(len(S)) :
            verArr = [None]*len( S[verlist] )
            for verInx in range( len( S[verlist] ) ) :
                verArr[verInx] = v[ S[verlist][verInx] ]
            segments[verlist] = verArr

        return segments
    
    # ---------------------------------------------------------------------+
    # Note:                                                                |
    # the following five _c...._vert methods assign color and values for   |
    # each vertex which are only used for FutDev export of vertex colors   |
    # (segment center colors are used for the matplotlib rendering)        |
    # ---------------------------------------------------------------------+ 

    def _map_cmap_from_sequence_vert(self,cmap) :
        i = np.linspace(0.0,1.0, len(self.vertexCoor) )    
        norm = colors.Normalize(vmin=0.0,vmax=1.0)

        color = cmap(norm(i))
        self.vertexColor = color
        self.vertexValues = list( range(len(self.vertexCoor)) )
        self._bounds['vertvlim'] = [ 0.0, 1.0 ]
        return

    def _map_cmap_from_direction_vert( self, edgecolors, edgevalues, vBounds ) :
        # since the direction is dependent on two vertices,
        # the segment colors are used for the vertexColors.
        # The last vertex color is assigned the last segment
        # segment color as an approximation.
        numEdges = len(edgecolors)

        lastColor = edgecolors[numEdges-1]
        vColor = copy.copy(edgecolors)
        vColor = np.append(vColor,[lastColor],axis=0)

        lastValue = edgevalues[numEdges-1]
        vValue = copy.copy(edgevalues)
        vValue = np.append(vValue,[lastValue],axis=0)
        
        self.vertexColor = vColor
        self.vertexValues = vValue
        self._bounds['vertvlim'] = vBounds
        return

    def _map_cmap_from_datagrid_vert(self,g, delta, dmin, cmap) :

        x,y,_ = np.transpose(self.vertexCoor)

        # normalize xy plane -1 to 1
        x = np.interp( x, (x.min(), x.max()), (-1.0, 1.0) )
        y = np.interp( y, (y.min(), y.max()), (-1.0, 1.0) )

        trial = [ g( x[i],y[i])[0] for i in range(len(x)) ]

        v = np.array(trial)*delta + dmin
        vmin, vmax = v.min(), v.max()
        norm = colors.Normalize(vmin=v.min(),vmax=v.max())

        color = cmap(norm(v))
        self.vertexColor = color
        self.vertexValues = v
        self._bounds['vertvlim'] = [ vmin, vmax ] 
        return 

    def _map_cmap_from_op_vert(self, operation, cmap) :
        xyz = np.transpose(self.vertexCoor)
        v = np.array(operation(xyz))
        vmin, vmax = v.min(), v.max()
        norm = colors.Normalize(vmin=v.min(),vmax=v.max())
        
        color = cmap(norm(v))
        self.vertexColor = color
        self.vertexValues = v
        self._bounds['vertvlim'] = [ vmin, vmax ] 
        return 

    def _map_color_from_op_vert(self, operation, rgb) :
        xyz = np.transpose(self.vertexCoor)
        colors = np.array(operation(xyz))
        colors = np.transpose(colors)
        colors = np.clip(colors,0,1)
        if rgb : RGB = colors
        else :   RGB = cm.colors.hsv_to_rgb(colors)
        if RGB.shape[1] == 3 :
            ones = np.ones(RGB.shape[0])[:,np.newaxis]
            RGB = np.concatenate((RGB,ones),axis=1)
        
        self.vertexColor = RGB
        self.vertexValues = None
        self._bounds['vertvlim'] = [ 0.0, 1.0 ] 
        return

    # ---------------------------------------------------------------------

    @property
    def _ledgecolors(self) :
        return self._edgecolors

    @_ledgecolors.setter
    def _ledgecolors(self,val) :
        self._edgecolors = val
        return

    @property
    def vlim(self) :
        '''Range of scalar values which may be associated with color'''
        return self._bounds['vlim']

    @vlim.setter
    def vlim(self,val) :
        self._bounds['vlim'] = val
        return

    @property
    def name(self) :
        '''Descriptive identifier for the line geometry.'''
        if self._geomName is None : return ''
        return self._geomName

    @name.setter
    def name(self, val) :
        self._geomName = val
        self.set_label( '' if val is None else val)
        return

    @property
    def cname(self) :
        '''Descriptive identifier for values indicated by color.'''
        if self.valuesName is None : return ''
        return self.valuesName

    @cname.setter
    def cname(self, val) :
        self.valuesName = val
        return

    @property
    def vertices(self) :
        """A 3 x N array of line vertices."""
        v = self.vertexCoor
        m = np.transpose(v)
        return m

    @property
    def segmentcenters(self) :
        """A 3 x N array of line segment centers."""  
        segments = self._segments3d
        totSegs = []
        for line in segments :
            npline = np.array(line)
            for segIx in range(npline.shape[0]-1) :
                center = np.add(npline[segIx],npline[segIx+1])/2
                totSegs.append(center)
        m = np.transpose(np.array(totSegs))
        return m

    @property
    def segmentdirections(self) :
        """A N x 3 array of line segment direction vectors.""" 
        # -----------------------------------------
        segments = self._segments3d
        if not self.isShredded :
            segments = self._shred()
        # -----------------------------------------
        vfT = np.transpose(segments, (1,0,2) )
        direction = vfT[1] - vfT[0]
        return direction

    @property
    def segmentcolors(self) :
        """A N x 4 array of segment colors."""
        return self._ledgecolors  # update for v_1.2.0

    @property
    def bounds(self) :
        """
        Dictionary of line geometric and value ranges.
        
        Each dictionary value is a 2 float array of minimum
        and maximum values of the line.  Keys are:

        'xlim' : x-coordinate

        'ylim' : y-coordinate

        'zlim' : z-coordinate

        'r_xy' : radial distance from the z axis

        'rorg' : radial distance from the origin

        'vlim' : value functional assignments.
            
        Values are assigned from the geometry and color mapping methods,
        including cllipping.

        """

        return self._bounds
    
    @property
    def cBar_ScalarMappable(self) :
        """matplotlib.cm.ScalarMappable object for line values. """

        return self.stcBar_ScalarMappable()

    def stcBar_ScalarMappable(self, stripe=None, bgcolor='k') :
        """
        ScalarMappable object with striped color values.
        Primarily useful for display of contour colorbars. 
        
        Parameters
        ----------
        stripe : integer, optional, default: None
            Number of color stripes.  If None, no
            stripes are formed.
        
        bgcolor : color, optional, default: 'black'
            Color between the stripes.
        
        Returns
        -------
        self : Matplotlib.cm.ScalarMappable object
        
        """
        # .......................................
        def stripeCmap( cmap, numb, background='k', cname=None) :
            if isinstance(cmap,str) : 
                cmap = cm.get_cmap(cmap)
            numbSegs = 256
            maxIndx = numbSegs-1
            bgrnd = cm.colors.to_rgba_array(background)
            cList = np.tile( bgrnd, (numbSegs,1))
            fIndex = np.linspace(0.0,1.0,numb)
            for i in fIndex :
                index = int(i*(maxIndx))
                cList[index] = cmap(i)
                if index < maxIndx : cList[index+1] = cmap(i)
                else : cList[index-2] = cmap(i)
                if index > 0 : cList[index-1] = cmap(i)
                else : cList[index+2] = cmap(i)
            return colors.ListedColormap(cList)
        # .......................................
        vmin, vmax = self._bounds['vlim']
        norm = colors.Normalize(vmin=vmin,vmax=vmax)
        objMap = self.get_cmap()
        if stripe is not None:
            objMap = stripeCmap( objMap, stripe, bgcolor) 
        sm = cm.ScalarMappable(cmap=objMap, norm=norm)
        sm.set_array([])
        return sm        

    def shred(self, rez=0) :
        """
        Recursively halve each line segment into separate lines.
        
        Parameters
        ----------
        rez : integer, optional, default: 0
            Number of recursive subdivisions of a line segment.
            Values range from 0 to 7.

        Returns
        -------
        self : line object
        
        """
        # .............................................................
        def midSegmentVertices(numDiv,tailCoor, headCoor) :
            # return added subdivision vertices, in order from
            # tail to head.
            minLen = 1.0/numDiv
            divArr = np.linspace(minLen,1.0-minLen, numDiv-1)
            delta = np.array(headCoor) - np.array(tailCoor)
            vSet = [ tailCoor + np.multiply(delta,mult) for mult in divArr  ]
            vSet = np.array(vSet)
            return vSet
        # .............................................................
        
        if not isinstance(rez, int) :
            raise ValueError('Incorrect rez type, must ba of type int')
        if (rez<0) or (rez>_MAXREZ) :
            raise ValueError('Incorrest rez of {}, must be an int, 0 <= rez <= {}'.format(rez,_MAXREZ))

        numDiv = int(2**rez)
        # FutDev: following seems to work in cases checked but still not 100% confident.
        if not self.isShredded :
            colors = self._shred_color()
            self._shred()
            self._ledgecolors = colors  # update for v_1.2.0
            # DEVNOTE: colors now refer to segment colors, not line colors.
        if rez==0: return self
        else :
           colors = self._shred_segm_color(numDiv) 

        numSeg = 0
        for line in self.lineIndices : numSeg += len(line)-1
        
        new_vertexCoor = np.copy(self.vertexCoor)
        new_lineIndices = []
        new_start_vert_index = len(self.vertexCoor)
        for line in self.lineIndices :
            single_lineIndex = []
            numLineSegs = len(line)-1
            lastLineVertIndex = numLineSegs
            for i in range(numLineSegs ) :
                tail_vert = self.vertexCoor[ line[i] ]   
                head_vert = self.vertexCoor[ line[i+1] ]
                single_lineIndex.append( line[i] )
                new_verts = midSegmentVertices(numDiv,tail_vert, head_vert)
                new_vertexCoor = np.concatenate( (new_vertexCoor,new_verts),axis=0)
                next_start_index = new_start_vert_index + len(new_verts)
                temp = [*range(new_start_vert_index,next_start_index)]
                single_lineIndex.extend( temp )
                new_start_vert_index = next_start_index
                if (i+1) == lastLineVertIndex : single_lineIndex.append( line[i+1] )
            new_lineIndices.append(single_lineIndex)

        self.vertexCoor = new_vertexCoor
        self.lineIndices = new_lineIndices
        self._shred()
        self._ledgecolors = colors  # update for v_1.2.0
        return self
    
    def map_cmap_from_sequence( self, cmap=None, sop=None, cname=None ) :
        """
        Sequential color assignment of each line segment.

        Parameters
        ----------
        cmap : str or Colormap, optional, default: 'viridis'
            A Colormap instance or registered colormap name.
        
        sop : function object, optional, default: None
            Function that takes one argument,
            a N Numpy array of scalar values ranging 0 to 1.
            The function returns a Numpy array of scalar values.

        cname : str, optional, default: function object name
                if not a lambda function, otherwise ''.
                Identifier for the color values.
        
        Returns
        -------
        self : line object

        """

        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cm_colorMap = self.get_cmap()

        self.valuesName = 'Sequence'   
        # -----------------------------------------
        segments = self._segments3d
        if not self.isShredded :
            segments = self._shred()
        # -----------------------------------------
        i = np.linspace(0.0,1.0, len(segments) )
        v = i
        vmin, vmax = 0.0, 1.0
        if sop is not None :
            temp = _getFunctionName(sop,cname)
            if temp is not None : self.valuesName = temp
            v = np.array(sop(i))
            vmin, vmax = v.min(), v.max()

        norm = colors.Normalize(vmin=vmin,vmax=vmax)
        self._ledgecolors = cm_colorMap(norm(v))  # update for v_1.2.0
        RGB = cm_colorMap(norm(v))
        self.set_color(RGB)

        self.segValues = v
        self._bounds['vlim'] = [ vmin, vmax ] 
        # .............................................
        self._map_cmap_from_sequence_vert(cm_colorMap)
        # .............................................
        return self   

    def map_cmap_from_direction( self, direction=None, cmap=None, isAbs=False, cname=None  ) :
        """
        Assignment of line color based on line segment direction.

        The dot product of line segment directions with the direction paramenter
        is used to assign line segment colors from a colormap.

        Parameters
        ----------
        direction : array-like, optional, default: (1,0,1)
            A xyz vector. 

        cmap : str or Colormap, optional, default: 'viridis'
            A Colormap instance or registered colormap name.
        
        isAbs : boolean, optional, default: False
                If True, the absolute value of the dot product
                will control the color.

        cname : str, optional, default is None
                Identifier for the color values.

        Returns
        -------
        self : line object

        """

        if direction is None : direction = _ILLUM

        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cm_colorMap = self.get_cmap()


        self.valuesName = cname
        # -----------------------------------------
        segments = self._segments3d
        if not self.isShredded :
            segments = self._shred()
        # -----------------------------------------
        segDir = self._get_segment_directions(segments)
        unitDirection = np.divide( direction, np.linalg.norm(direction) )
        v = np.dot(segDir,unitDirection)
        if isAbs : v = -np.abs(v)
        vmin, vmax = v.min(), v.max()
        norm = colors.Normalize(vmin=v.min(),vmax=v.max())

        self._ledgecolors = cm_colorMap(norm(v))  # update for v_1.2.0
        RGB = cm_colorMap(norm(v))
        self.set_color(RGB)

        self.segValues = v
        self._bounds['vlim'] = [ vmin, vmax ] 
        # .............................................
        self._map_cmap_from_direction_vert( self._ledgecolors, self.segValues, self._bounds['vlim'] )  # update for v_1.2.0
        # .............................................
        return self

    def map_cmap_from_datagrid(self, datagrid, cmap=None, cname=None ) :
        """
        Line color assignment using a 2D datagrid.

        Datagrid values are normalized in the range 0 to 1.

        Parameters
        ----------
        datagrid : 2D float array

        cmap : str or Colormap, optional
            A Colormap instance or registered colormap name.
            If not assigned, the surface Colormap is used.
            The colormap maps the datagrid values to colors.

        cname : string, optional, default: None.
            Descriptive identifier for values indicated by color.

        Returns
        -------
        self : line object

        """

        data = datagrid
        dmax = np.amax(datagrid)
        dmin = np.amin(datagrid)
        delta = dmax-dmin
        data = (datagrid - dmin)/delta
        xd = np.linspace(-1, 1, data.shape[0] )
        yd = np.linspace(-1, 1, data.shape[1] )
        g =  interpolate.interp2d(yd, xd, data, kind='cubic')

        # -----------------------------------------
        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cm_colorMap = self.get_cmap()

        self.valuesName = cname

        # -----------------------------------------
        segments = self._segments3d
        if not self.isShredded :
            segments = self._shred()
        # -----------------------------------------
        segCoor = self._get_segment_centers(segments)
        x,y,_ = np.transpose(segCoor)

        # normalize xy plane -1 to 1
        x = np.interp( x, (x.min(), x.max()), (-1.0, 1.0) )
        y = np.interp( y, (y.min(), y.max()), (-1.0, 1.0) )

        trial = [ g( x[i],y[i])[0] for i in range(len(x)) ]

        v = np.array(trial)*delta + dmin
        vmin, vmax = v.min(), v.max()
        norm = colors.Normalize(vmin=v.min(),vmax=v.max())

        self._ledgecolors = cm_colorMap(norm(v))  # update for v_1.2.0
        RGB = cm_colorMap(norm(v))
        self.set_color(RGB)

        self.segValues = v
        self._bounds['vlim'] = [ vmin, vmax ] 

        # .............................................
        self._map_cmap_from_datagrid_vert(g, delta, dmin, cm_colorMap)
        # .............................................

        return self

    def map_cmap_from_op( self, operation, cmap=None, cname=None ) :
        """
        Functional assignment of a color from a color map.

        Line segment center coordinates are used to calculate a scalar which
        is then used to assign segment colors from a colormap.

        Parameters
        ----------
        operation : function object
            Function that takes one argument,
            a 3xN Numpy array of native coordinates.
            The function returns a Numpy array of scalar values.

        cmap : str or Colormap, optional, default: 'viridis'
            A Colormap instance or registered colormap name.
        
        cname : str, optional, default: function object name
                if not a lambda function, otherwise ''.
                Identifier for the color values.

        Returns
        -------
        self : line object

        """       
        cname = _getFunctionName(operation,cname)
        if cname is not None : self.valuesName = cname

        if cmap is not None :
            if isinstance(cmap,str) : 
                self.set_cmap(cm.get_cmap(cmap))
            else :
                self.set_cmap(cmap)
        cm_colorMap = self.get_cmap()

        # -----------------------------------------
        segments = self._segments3d
        if not self.isShredded :
            segments = self._shred()
        # -----------------------------------------
        segCoor = self._get_segment_centers(segments)
        # .....
        xyz = np.transpose(segCoor)
        v = np.array(operation(xyz))
        vmin, vmax = v.min(), v.max()
        norm = colors.Normalize(vmin=v.min(),vmax=v.max())
        
        self._ledgecolors = cm_colorMap(norm(v))  # update for v_1.2.0
        RGB = cm_colorMap(norm(v))
        self.set_color(RGB)

        self.segValues = v
        self._bounds['vlim'] = [ vmin, vmax ] 
        # ========================================================================
        self._map_cmap_from_op_vert(operation, cm_colorMap)
        # ========================================================================
        return self    

    def map_color_from_op( self, operation, rgb=True, cname=None ) :
        """
        Assignment of line color from a function.

        Line segment colors are assigned from a function of line  
        segment center coordinates.

        Parameters
        ----------
        operation : function object
            Function that takes one argument,
            a 3xN Numpy array of native coordinates.
            The function returns a 3xN color value.

        rgb : bool {True, False}, optional, default: True
            By default, RGB color values are returned by the
            operation function.  If set False, the operation
            returns HSV color values.
        
        cname : str, optional, default: function object name
                if not a lambda function, otherwise ''.
                Identifier for the color values.
                
        Returns
        -------
        self : line object

        """       
        cname = _getFunctionName(operation,cname)
        if cname is not None : self.valuesName = cname
        # -----------------------------------------
        segments = self._segments3d
        if not self.isShredded :
            segments =self._shred()
        # -----------------------------------------
        segCoor = self._get_segment_centers(segments)

        xyz = np.transpose(segCoor)
        colors = np.array(operation(xyz))
        colors = np.transpose(colors)
        colors = np.clip(colors,0,1)
        if rgb : RGB = colors
        else :   RGB = cm.colors.hsv_to_rgb(colors)
        if RGB.shape[1] == 3 :
            ones = np.ones(RGB.shape[0])[:,np.newaxis]
            RGB = np.concatenate((RGB,ones),axis=1)
        
        self._ledgecolors = RGB  # update for v_1.2.0
        self.set_color(RGB)

        self.segValues = None
        self._bounds['vlim'] = [ 0.0, 1.0 ] 
        # ========================================================================
        self._map_color_from_op_vert( operation, rgb )
        # ========================================================================
        return self    

    def map_geom_from_op(self, operation, name=None) :
        """
        Functional transformation of line vertex coordinates.

        Parameters
        ----------
        operation : function object
            Coordinate mapping function that takes one
            argument, a 3xN Numpy array of native coordinates.
            The function returns a 3xN array of coordinates.

        name : string, optional, default: None or op function name.
            Descriptive identifier for the geometry.

        Returns
        -------
        self : line object

        """
        orgCoor = self.vertexCoor
        xyz = np.transpose(orgCoor)
        XYZ = np.array(operation(xyz))
        v = np.transpose(XYZ)
        self.vertexCoor = v
        segments = self._get_segments( self.vertexCoor, self.segmIndices  )
        self.set_segments(segments)
        self._set_geometric_bounds()
        if name is None : name = self._geomName
        name = _getFunctionName(operation,name)
        if name is not None : self._geomName = name
        return self

    def map_to_plane(self, dist, **kargs) :
        """
        Project line onto a mapping surface.

        Parameters
        ----------
        dist : number, optional, default : 0.0 
            distance from the origin to the surface,
            along the direction vector.

        direction : array-like, optional, default: [0,0,1]
            A xyz vector normal to the intersection planes
            for planar mapping or axial direction of 
            a cylinder for cylindrical mapping.

        coor : integer or string indicating the type of mapping surface:
            0, p, P, xyz,planar                 - planar (default)
            1, c, C, cylinder,pplar,cylindrical - cylinder
            2, s, S, sphere,spherical           - sphere
            3, x, X                             - y-z plane
            4, y, Y                             - x-z plane
            5, z, Z                             - x-y plane

        Returns
        -------
        self : line object

        """
        dirDft = { 'direction': [0,0,1.0] }
        kargDft = { **_COOR_KWARGS, **dirDft }
        KW = KWprocessor( kargDft, 'map_to_plane', **kargs)
        coor = KW.getVal('coor')
        direction = KW.getVal('direction')

        if coor==3 : coor, direction = 0, [ 1,0,0 ]
        if coor==4 : coor, direction = 0, [ 0,1,0 ]
        if coor==5 : coor, direction = 0, [ 0,0,1 ]

        # .......................................................
        def project_plane(xyz) :
            verts = xyz.T
            unitDirection = np.divide( direction, np.linalg.norm(direction) )
            dprod = np.dot( verts, unitDirection )
            vdot = np.array([dprod]).T
            dUnit = np.tile(unitDirection, (len(verts),1) )
            vUnit = vdot*dUnit
            dUnit = dist*dUnit
            projVerts = verts + (dUnit - vUnit)
            return projVerts.T
        def project_cyl(xyz) :
            verts = xyz.T
            unitDirection = np.divide( direction, np.linalg.norm(direction) )
            dmag = np.dot(verts,unitDirection)
            uVect = np.tile(unitDirection,(len(dmag),1) )
            dVect = np.array([dmag]).T
            A = uVect*dVect
            zeta = verts - A
            zetanorm = np.linalg.norm(zeta,axis=1)
            znorm = np.tile(zetanorm, (3,1))
            unitZeta = np.divide( zeta, znorm.T) 
            # above may report nan warning, but following corrects these values.
            unitZeta[ zetanorm <= 0.0, : ] = [0.0,0.0,0.0]
            projVerts = A + dist*unitZeta
            return projVerts.T
        def project_sph(xyz) :
            verts = xyz.T
            dprod = np.linalg.norm(verts,axis=1)
            tnorm = np.tile(dprod, (3,1) )
            unitVerts = np.divide( verts, tnorm.T )
            projVerts = dist * unitVerts
            return projVerts.T
        # .......................................................
        project = project_plane
        if coor == 1 : project = project_cyl
        if coor == 2 : project = project_sph

        return self.map_geom_from_op( project )

    def clip(self, operation) :
        """
        Remove segments from the line.
        NOTE: all operations are in xyz coordinates.

        Parameters
        ----------
        operation : function object
            Function that takes one argument, a 3xN Numpy array of coordinates.
            The function returns N array of bool { True, False } indicating if
            the segment, at the segment center coordinate, is to be 
            retained (True).  Otherwise, the segment is removed from the
            line (False).

        Returns
        -------
        self : line object

        """
        # -----------------------------------------
        colors = self._ledgecolors  # update for v_1.2.0
        segments = self._segments3d

        if self.restoredClip is None :
            self.restoredClip = { 
                'indices' : self.segmIndices, 
                'colors' : colors,
                'isShredded' : self.isShredded }

        if not self.isShredded :
            # .... note:colors modified, not created
            colors = self._shred_color()
            segments = self._shred()
            self._ledgecolors = colors  # update for v_1.2.0
        else :
            lenCol = len(self._ledgecolors)  # update for v_1.2.0
            lenSeg = len(segments)
            if lenCol != lenSeg : colors = self._shred_color()
        # -----------------------------------------
        segCoor = self._get_segment_centers(segments)
        xyz = np.transpose(segCoor)
        shouldKeep = np.array(operation(xyz))

        if np.any(shouldKeep) is None :
            warnings.warn('WARNING: Clipping resulted in no lines found, line return unclipped.')
            return self

        clipIndices = np.array(self.segmIndices)[shouldKeep]
        clipColors = colors[shouldKeep]
        
        clipColors = np.array(clipColors)
        clipCoor = self.vertexCoor[  np.array(clipIndices) ]   # the lineset is shredded so can do this.
        self.set_segments(clipCoor)
        self.set_color(clipColors)

        # FutDev: 'correct' definition should use nparrays for segmIndices not lists 
        # self.segmIndices = clipIndices  ... bug fix,  just a bandaid.
        self.segmIndices = clipIndices.tolist()

        return self

    def clip_plane(self,dist,**kargs) :
        """
        Remove segments from the line based on a clip surface.

        Parameters
        ----------
        dist : number, optional, default : 0.0 
            Distance from the origin to the clip surface,
            along the direction vector.

        direction : array-like, optional, default: [0,0,1]
            A xyz vector normal to the intersection plane
            for a planar clip surface or axial direction of 
            a cylinder for a cylindrical clip surface.

        coor : integer or string indicating the type of clip surface:
            0, p, P, xyz,planar                 - planar (default)
            1, c, C, cylinder,pplar,cylindrical - cylinder
            2, s, S, sphere,spherical           - sphere
            3, x, X                             - y-z plane
            4, y, Y                             - x-z plane
            5, z, Z                             - x-y plane

        Returns
        -------
        self : line object

        """
         # ...........................................................
        def clip_op(xyz,dist,coor,direction) :
            abc = np.transpose(xyz)
            dotprod = _coor_dotprod(direction,coor,abc)
            shouldclip = np.less(dotprod,np.full(len(dotprod),dist))
            return shouldclip
        # ...........................................................
        dirDft = { 'direction': [0,0,1.0] }
        kargDft = { **_COOR_KWARGS, **dirDft }
        KW = KWprocessor( kargDft, 'clip_plane', **kargs)
        coor = KW.getVal('coor')
        direction = KW.getVal('direction')

        if coor==3 : coor, direction = 0, [ 1,0,0 ]
        if coor==4 : coor, direction = 0, [ 0,1,0 ]
        if coor==5 : coor, direction = 0, [ 0,0,1 ]
        
        return self.clip(lambda xyz : clip_op(xyz,dist,coor,direction) )

    def _restore_from_clip(self) :
        """
        Reset segment and colors prior to any clip operation.
        (in development)
        """
        
        if self.restoredClip is None :
            warnings.warn('There are no clipped lines to restore.')
            return self

        self.isShredded = self.restoredClip['isShredded']
        self.set_color( self.restoredClip['colors'] )

        self.segmIndices = self.restoredClip['indices']
        segmIndices = self.segmIndices

        segments = self._get_segments( self.vertexCoor, segmIndices  )
        self.set_segments(segments)    
        self.restoredClip = None
        return self

    def set_line_alpha(self, alpha, constant=False) :
        """
        Adjust the alpha value of the line segment colors.

        Parameters
        ----------
        alpha : scalar
              Alpha is in the domain 0 to 1.

        constant : bool { True, False }, optional, False
            If False, segment color values are multiplied by
            alpha.  If True, all segment colors alpha channels
            are assigned to alpha.

        Returns
        -------
        self : line object

        """
        
        if (alpha<0.0) or (alpha>1.0) :
            raise ValueError('line alpha must be between 0 and 1, found {}'.format(alpha))

        orig_colors = self._ledgecolors  # update for v_1.2.0
        if constant :
            np.put_along_axis(orig_colors, np.array([[3]]), alpha, axis=1)
            colorMap = orig_colors
        else:
            colorMap = np.multiply(orig_colors,np.array([1.0,1.0,1.0,alpha]))
        
        self._ledgecolors = colorMap  # update for v_1.2.0

        return self

    def shade(self, depth=0, direction=None, contrast=None) :
        """
        Reduce line HSV color Value based on line segment direction.
        
        The dot product of line segment with the illumination direction 
        is used to adjust face HSV color value.

        Parameters
        ----------
        depth : scalar, optional, default: 0
            Minimum color value of shaded line segments.
            Depth value ranges from 0 to 1.

        direction : array-like, optional, default: (1,0,1)
            A xyz vector pointing to the illumination source. 

        contrast : scalar, optional, default: 1
            Shading contrast adjustment from low to high with a value
            of 1 for linear variations with the line segment direction.
            Contrast value ranges from 0.1 to 3.     
        
        Returns
        -------
        self : line object

        """

        if direction is None : direction = _ILLUM
        if np.any( (depth<0) | (depth>1) ) :
            raise ValueError('depth values, {}, must be between 0 and 1'.format(depth))
        if contrast is not None:
            if (contrast<0.1 or contrast>3) :
                raise ValueError('contrast must be between 0.1 and 3. , found {}'.format(contrast))


        # -----------------------------------------
        colors = self._ledgecolors  # update for v_1.2.0
        segments = self._segments3d
        if not self.isShredded :
            # .... note:colors modified, not created
            colors = self._shred_color()
            segments = self._shred()
            self._ledgecolors = colors  # update for v_1.2.0
        else :
            lenCol = len(self._ledgecolors)  # update for v_1.2.0
            lenSeg = len(segments)
            if lenCol != lenSeg : colors = self._shred_color()
        # -----------------------------------------
        segDir = self._get_segment_directions(segments)
        unitDirection = np.divide( direction, np.linalg.norm(direction) )
        dtprod = np.dot(segDir,unitDirection)
        # for Lines, effective shade is normal to line
        d = 1- np.abs(dtprod)
        if contrast is not None :
            d =  0.5*( 1 + np.sign(dtprod)*np.power( np.abs(dtprod) , 1/contrast) )      
        orig_colors = colors
        alphas = orig_colors[:,3:4]
        fc_less_alpha = orig_colors[:,:3]
        hsv_vals = cm.colors.rgb_to_hsv(fc_less_alpha)
        # adjust color value with multiplier.
        cvm = ((1-depth)*d + depth*np.ones(len(d)))[:, np.newaxis]
        hsv_vals[:,2] = (cvm*hsv_vals[:,2:3])[:,0]
        rgb_vals =  cm.colors.hsv_to_rgb(hsv_vals)
        shade_colors = np.concatenate((rgb_vals, alphas), axis=1)
        self._ledgecolors = shade_colors  # update for v_1.2.0
        self.set_color(shade_colors)

        return self

    def fade(self,depth=0,elev=None,azim=None,ax=None) :
        """
        Reduce line opacity based on line segment position relative
        to the view orientation.
        
        Parameters
        ----------
        depth : scalar, optional, default: 0
            Minimum opacity to 1 for linear line opacity from back
            to front line segment position.
            Depth value ranges from 0 to 1.

        elev : scalar, optional, default: 30
            Elevation of the axes view.

        azim : scalar, optional, default: -60
            Azimuth of the axes view. 
        
        ax: Matplotlib 3D axes.
            If not None, elev and azim are taken from the ax.
        
        Returns
        -------
        self : line object

        """

        if ax is None :
            if elev is None: elev = _DFT_VIEW[0]
            if azim is None: azim = _DFT_VIEW[1]
        else :
            elev = ax.elev
            azim = ax.azim
        direction = elev_azim_2vector(elev, azim)
        unitDirection = np.divide( direction, np.linalg.norm(direction) )

        if np.any( (depth<0) or (depth>1) ) :
            raise ValueError('depth values, {}, must be between 0 and 1'.format(depth))

        # -----------------------------------------
        #colors = self._edgecolors
        colors = self.get_edgecolor()  # update for v_1.2.0
        segments = self._segments3d
        if not self.isShredded :
            # .... note:colors modified, not created
            colors = self._shred_color()
            segments = self._shred()
            #self._edgecolors = colors
            self.set_edgecolor(colors)  # update for v_1.2.0
        else :
            #lenCol = len(self._edgecolors)
            lenCol = len(colors)  # update for v_1.2.0
            lenSeg = len(segments)
            if lenCol != lenSeg : colors = self._shred_color()
        # -----------------------------------------
        segCenters = self._get_segment_centers(segments)
        if len(segCenters) == 1 : return self # fade not available for single segment line 
        dtprod = np.dot(segCenters,unitDirection)
        dprange = np.amax(dtprod)-np.amin(dtprod)
        norm_dp = (dtprod - np.amin(dtprod))/dprange
        fadex = (1-depth)*norm_dp + depth

        orig_colors = colors
        fc_less_alpha = orig_colors[:,:3]
        fadex = fadex*orig_colors[:,3]  # update for v_1.2.0 
        faded_colors = np.concatenate((fc_less_alpha, fadex[:,np.newaxis] ), axis=1)
        #self._edgecolors = faded_colors   
        self.set_edgecolor(faded_colors)  # update for v_1.2.0   
        self.set_color(faded_colors)  # update for v_1.2.0   

        return self

    def transform(self, rotate=None, scale=1.0, translate=[0,0,0] )  :
        """
        Linear transformation of the line object.

        Parameters
        ----------
        rotate :  array-like, optional, default: None
            A 3 by 3 rotation matrix.

        scale : 3D array or scalar, optional, default: 1.0
            Multipliers of the object coordinates, in
            xyz coordinates. If a single scalar, three
            multipliers are the same value.
        
        translate : array-like, optional, default: [0,0,0]
            A xyz translation vector of object origin. 

        Returns
        -------
        self : line object

        """
        if rotate is None : rotate = np.identity(3)
        if isinstance(scale,(int,float)) :
            scale =  np.full(3,scale).astype(float)
        #.....................................    
        def isOrth(M) :
            M = np.array(M)
            tol = 0.00001
            prod = np.dot(M,M.T)
            np.fill_diagonal(prod,0.0)
            return not (np.abs(prod)>tol).any()
        #.....................................    
        trans = lambda V,S,T : np.add(  np.dot( np.multiply(V,S),rotate) ,T)
        #.....................................    
        if not isOrth(rotate) :
            warnings.warn('Transform rotate matrix is NOT Orthoginal.')
        
        v = trans( self.vertexCoor, scale, translate  )
        segs = self._get_segments( v, self.segmIndices)

        self.vertexCoor = v
        self.set_segments(segs)
        self.scale = scale
        self.translation = translate
        self.rotation = rotate

        _set_geomBounds(self.vertexCoor,self._bounds)

        return self

    def get_transformAxis(self, **kargs) :
        """
        Line3DCollection of the 'last' transformed line coordinate axis.
        
        Parameters
        ----------
        lenmult : scalar or 3D array, optional, default: 1
            Scalar multiplier of the three coordinate axis.

        width : scalar, optional, default: 1.5
            Line width of the coordinate axis.   

        color : a Color or 3D array of Colors, optional, default: ['r','g','b']

        negaxis : boolean, optional, default: False
            If True, include the negative axis, otherwise axes start at origin.

        Returns
        -------
        Line3DCollection object

        """
        kargDft = { 'lenmult':1.0, 'width':1.5, 'color':['r','g','b'], 'negaxis':False }
        KW = KWprocessor(kargDft,'get_transformAxis',**kargs)
        lenmult = KW.getVal('lenmult')
        width = KW.getVal('width')
        color = KW.getVal('color')
        negaxis = KW.getVal('negaxis')

        selfProp = {
            'scale' : self.scale,
            'rotation' : self.rotation,
            'translation' : self.translation
        }
        return _transformAxis(selfProp,lenmult, width, color, negaxis)

    def _rectangulateFilledSurFace(self,rez,baseVcoor,baseFaceVertexIndices, midVectFunc) :
        vCoor = baseVcoor.copy()
        fvIndices = []
        evIndices = []
        #.....................................    
        def getVerticesLine(degree, leftIndex, rightIndex) :
            def recurs_edgeCoor(degree, indexOrigin, isRight, leftIndex, rightIndex, Eindex) :
                m = degree - 1
                delta = int(np.power(2,m))
                if isRight : delta = -delta
                i = indexOrigin - delta
                mid_coor = midVectFunc( vCoor[leftIndex], vCoor[rightIndex] )
                currentCoorIndex = len(vCoor)
                vCoor.append(mid_coor)
                Eindex[i] = currentCoorIndex
                if m <= 0 :
                    evIndices.append( [leftIndex, currentCoorIndex] )                    
                    evIndices.append( [currentCoorIndex, rightIndex] )                    
                    return
                recurs_edgeCoor(m,i,False,leftIndex,currentCoorIndex, Eindex)
                recurs_edgeCoor(m,i,True,currentCoorIndex,rightIndex, Eindex)
                return
            #.....................................    
            iOrigin = int(np.power(2,degree))
            Eindex = [None]*(iOrigin+1)
            Eindex[0] = leftIndex
            Eindex[iOrigin] = rightIndex
            if degree==0 :
                evIndices.append( [leftIndex, rightIndex] )
                return Eindex
            recurs_edgeCoor(degree, iOrigin, False, leftIndex, rightIndex, Eindex)
            return Eindex
        #.....................................    
        def recurs_faceIndices(A,C,rez):
            # similar to inner function of triangulateBase
            if rez==0 :
                abc = [ A[0], A[1], C[0], C[1] ]               
                fvIndices.append( abc )
                return
            # -----------------------------------------------------
            J = int( (len(A)-1)/2 )
            K = J + 1
            # sudivide the two inner quadraterals .................
            recurs_faceIndices( A[:K],  C[J:], rez-1)
            recurs_faceIndices( C[:K],  A[J:], rez-1)
            return
        #.....................................
        i_prev_2 = -1  # check for first face only.
        for face in baseFaceVertexIndices :
            if i_prev_2 == face[3] :
                B = D[::-1]
            else :
                B = getVerticesLine(rez, face[3],face[0] )
            i_prev_2 =  face[2]
            D = getVerticesLine(rez, face[1],face[2] )
            recurs_faceIndices(B,D,rez)

        vertexCoor = np.array(vCoor) 
        indexObj = { 'face': np.array(fvIndices) , 'edge': np.array(evIndices) }
        return indexObj ,vertexCoor

    def _get_filled_surface(self, addedverts, lrez=0, sclr=None, alpha=None, name=None) :
        # Dual purpose algorithm for lines projected to a plane (True) or
        # Parametric lines projected to another Parametric line.

        #.....................................
        def prepSegs() :
            # devNote: even if lines are not shredded, each
            # line could possibly have different colors.
            # Hence, must construct surfaces from 'shredded'
            # lines without actually shredding the line.
            # The color array is then for each polygon.
            # The following is actually somewhat a duplication
            # of the _shred.
            if not self.isShredded :
                color = self._shred_color()
                linesegs = self.lineIndices
                indices = []
                for line in linesegs :
                    numSegs = len(line) - 1
                    for seg in range (numSegs) :
                        t = [ line[seg],line[seg+1] ]
                        indices.append(t)
                vertIx = np.array(indices)
            else :
                color = self._edgecolors  # update for v_1.2.0
                vertIx= np.array(self.segmIndices)
            return vertIx, color
        #.....................................
        
        topSegInx,color = prepSegs()

        stInx = len(self.vertexCoor)
        botSegInx = topSegInx + stInx
        topSegInx = np.flip(topSegInx,axis=1)
        fcIndx = np.concatenate( (botSegInx,topSegInx), axis=1)
           
        sverts = np.concatenate((self.vertexCoor,addedverts),axis=0)
        scolor = np.array(color)

        if lrez > 0 :
            sverts = list(sverts)
            indexObj ,vertexCoor = self._rectangulateFilledSurFace(lrez,sverts,fcIndx, self._midVectorFun)
            sverts = vertexCoor
            fcIndx = indexObj['face']
            scolor = np.repeat(scolor,2**lrez,axis=0)

        if sclr is not None : scolor = sclr

        if alpha is not None:
            surface = Surface3DCollection(sverts, fcIndx, name=name, color=scolor, alpha=alpha)
        else :
            surface = Surface3DCollection(sverts, fcIndx, name=name,color=scolor)
        return surface

    def get_filled_surface(self,**kargs) :
        """
        Surface3DCollection from a line projected to a surface.

        Parameters
        ----------
        dist : number, optional, default : 0.0 
            distance from the origin to the projected surface,
            along the direction vector.

        direction : array-like, optional, default: [0,0,1]
            A xyz vector normal to the intersection plane
            for planar projections or axial direction of 
            a cylinder for cylindrical projections.

        coor : integer or string indicating the type of projection surface:
            0, p, P, xyz,planar                 - planar (default)
            1, c, C, cylinder,pplar,cylindrical - cylinder
            2, s, S, sphere,spherical           - sphere
            3, x, X                             - y-z plane
            4, y, Y                             - x-z plane
            5, z, Z                             - x-y plane

        lrez : positive integer, optional, default : 0
            projection recursive subdivisions.

        color : Color, optional, default : line color

        alpha : scalar, optional, default : 1
            Face color alpha value, in the range 0 to 1. 
            
        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Returns
        -------
        Surface3DCollection object
       
        """
        extDft = { 'direction':[0,0,1], 'dist':0.0, 'lrez':0, 'alpha':None, 'color':None, 'name':None }
        kargDft = { **_COOR_KWARGS, **extDft }
        KW = KWprocessor(kargDft,'get_filled_surface',**kargs)
        drctn = KW.getVal('direction')
        axis = KW.getVal('coor')
        base = KW.getVal('dist')
        lrez = KW.getVal('lrez')
        alpha = KW.getVal('alpha')
        color = KW.getVal('color')
        name = KW.getVal('name')

        if base is None : base = 0.0
        if axis==3 : axis, drctn = 0, [ 1,0,0 ]
        if axis==4 : axis, drctn = 0, [ 0,1,0 ]
        if axis==5 : axis, drctn = 0, [ 0,0,1 ]

        selfcopy = copy.copy(self)
        mappedself = selfcopy.map_to_plane( base, direction=drctn, coor=axis)
        surf = self.get_surface_to_line(mappedself, lrez=lrez, alpha=alpha, color=color, name=name)
        return surf

    def get_surface_to_line(self, line, **kargs) :
        """
        Surface3DCollection from a line projected to another line.

        Parameters
        ----------
        line : line object with the same number of vertices
            as the calling object.

        lrez : positive integer, optional, default : 0
            projection recursive subdivisions.

        color : Color, optional, optional, default : line color

        alpha : scalar, optional, default : None (1)
            Face color alpha value, in the range 0 to 1. 
            
        name : string, optional, default: None.
            Descriptive identifier for the geometry.

        Returns
        -------
        Surface3DCollection object
       
        """

        kargDft = { 'lrez':0, 'alpha':None, 'color':None, 'name':None }
        KW = KWprocessor(kargDft,'get_surface_to_line',**kargs)
        lrez = KW.getVal('lrez')
        alpha = KW.getVal('alpha')
        name = KW.getVal('name')
        color= KW.getVal('color')

        len1, len2 = len(self.vertexCoor), len(line.vertexCoor)
        if len1 != len2 :
            raise ValueError('Mismaatched line rez for Parametric Lines: {} != {}'.format(len1,len2))

        addedverts = line.vertexCoor

        return self._get_filled_surface(addedverts, lrez=lrez, sclr=color, alpha=alpha, name=name)


class ParametricLine(ColorLine3DCollection) :
    
    def __init__( self, rez , operation=None, name=None, **kwargs ) :
        """
        Create single line from a parametric function with domain [0,1]

        Parameters
        ----------
        rez : integer, optional, default: 0 
            Number of recursive bisection of line segments.
            Rez values range from 0 to 7. If rez is less than zero,
            negative value will be the number of segments.

        operation : function object 
            Function that takes one argument, a Numpy float
            array. The function domain is [0.0, 1,0] and returns
            a 3 X N array of xyz coordinates.

        name : str, optional, default: function object name
                if not a lambda function, otherwise string ''.
                
        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.ColorLine3DCollection'.

        """
        # -----------------------------------------
        def default_line(numVerts) :
            # xyz axes lines [-1,1]
            lineVals = np.linspace(-1.0,1.0, int(numVerts) )
            zeros = np.zeros( int(numVerts) )
            Xaxis = np.asarray( [lineVals, zeros, zeros ] ).T
            Yaxis = np.asarray( [zeros, lineVals, zeros ] ).T
            Zaxis = np.asarray( [zeros, zeros, lineVals ] ).T
            verts = np.vstack( (Xaxis,Yaxis,Zaxis ))
            indices = [None]*3
            index = 0
            for i in range(0,3) :
                indices[i] = [None]*(numVerts)
                for j in range(numVerts) :
                    indices[i][j] = index
                    index += 1
            return verts, indices
        # -----------------------------------------
        if not isinstance(rez,int) :
            raise ValueError('ParametricLine: Incorrest rez of {}, must be an int'.format(rez))
        if (rez>_MAXPLREZ) :
            raise ValueError('ParametricLine:Incorrest rez of {}, must less than or equal to {}'.format(rez,_MAXPLREZ))
        if rez<0 :
            self.numSeg = -rez
            # FutDev: should check if > 2**(_MAXPLREZ+1)
        else :
            self.numSeg = int(2**(rez+1))
        numVerts = self.numSeg + 1
        if operation is None :
            verts, indices = default_line(numVerts)
        else :
            name = _getFunctionName(operation,name)
            lineVals = np.linspace(0.0,1.0, self.numSeg+1)
            verts = np.array(operation(lineVals)).T
            verts = verts.tolist()
            indices = [ [*range(len(verts))] ]
        super().__init__(verts, indices, name, **kwargs)

        self.rez = rez

    def _getSlice(self, g, numVerts, domain, xVal=None, yVal=None  ) :
        if xVal is None and yVal is None : return None, None
        numLines = 0
        if xVal is not None : numLines += 1
        if yVal is not None : numLines += 1           
        indices = [None]*( numLines )
        numb_X_lines = 0
        index = 0

        verts = []
        if xVal is not None :
            numb_X_lines = 1
            lineVertsPts = np.linspace(domain['ymin'], domain['ymax'], numVerts )
            x_lineVertsPts = np.array([xVal] )
            z_X_lines = g(x_lineVertsPts,lineVertsPts)

            for iY in range(numVerts) :
                 verts.append( [ xVal , lineVertsPts[iY], z_X_lines[iY,0] ])

            indices[0] = [None]*(numVerts)
            for j in range(numVerts) :
                indices[0][j] = index
                index += 1
        
        if yVal is not None :
            lineVertsPts = np.linspace(domain['xmin'], domain['xmax'], numVerts )
            y_lineVertsPts = np.array([yVal] )
            z_Y_lines = g(lineVertsPts, y_lineVertsPts)
            # the following, compensates for interp2d for datagrid
            # which drops the first dimension, see map_xySlice_from_op  ???
            z_Y_lines = np.array([g(lineVertsPts, y_lineVertsPts)])
            
            for iX in range(numVerts) :
                 verts.append( [ lineVertsPts[iX], yVal, z_Y_lines[0,iX] ])
            
            indices[numb_X_lines] = [None]*(numVerts)
            for j in range(numVerts) :
                indices[numb_X_lines][j] = index
                index += 1

        verts = np.array(verts)

        return verts, indices

    def _getSliceSet(self, g, numVerts, domain, xsamples, ysamples ) :
        
        if xsamples is None and ysamples is None : return None, None
        numb_X_lines, numb_Y_lines = 0,0
        if xsamples is not None : numb_X_lines = xsamples
        if ysamples is not None : numb_Y_lines = ysamples
        indices = [None]*( numb_X_lines + numb_Y_lines )
        index = 0

        verts = []
        if xsamples is not None :
            lineVertsPts = np.linspace(domain['ymin'], domain['ymax'], numVerts )
            x_lineVertsPts = np.linspace(domain['xmin'], domain['xmax'], numb_X_lines )
            if xsamples == 1 : x_lineVertsPts = np.asarray([0.0])
            z_X_lines = g(x_lineVertsPts,lineVertsPts)

            for iX in range(numb_X_lines) :
                for iY in range(numVerts) :
                    verts.append( [ x_lineVertsPts[iX] , lineVertsPts[iY], z_X_lines[iY,iX] ])

            for i in range(0,numb_X_lines) :
                indices[i] = [None]*(numVerts)
                for j in range(numVerts) :
                    indices[i][j] = index
                    index += 1

        if ysamples is not None :
            lineVertsPts = np.linspace(domain['xmin'], domain['xmax'], numVerts )
            y_lineVertsPts = np.linspace(domain['ymin'], domain['ymax'], numb_Y_lines )
            if ysamples == 1 : y_lineVertsPts = np.asarray([0.0])
            z_Y_lines = g(lineVertsPts,y_lineVertsPts)

            for iY in range(numb_Y_lines) :
                for iX in range(numVerts) :
                    verts.append( [  lineVertsPts[iX], y_lineVertsPts[iY] ,z_Y_lines[iY,iX] ])

            for i in range(numb_X_lines,numb_X_lines+numb_Y_lines) :
                indices[i] = [None]*(numVerts)
                for j in range(numVerts) :
                    indices[i][j] = index
                    index += 1
        
        verts = np.array(verts)

        return verts, indices

    def _map_xySlice(self, g, **kargs ) :       
        xplane, yplane, xsamples, ysamples = None,None,None,None
        if 'xplane' in kargs : xplane = kargs['xplane']
        if 'yplane' in kargs : yplane = kargs['yplane']
        if 'xset' in kargs : xsamples = kargs['xset']
        if 'yset' in kargs : ysamples = kargs['yset']

        planeNormal = None  # both x and y slice planes are defined.
        if (xplane is None) and (xsamples is None) :
            planeNormal = [0,1,0]  # must be in the x-z plane
        if (yplane is None) and (ysamples is None) :
            planeNormal = [1,0,0]  # must be in the y-z plane

        domain = {'xmin':-1, 'xmax':1, 'ymin':-1, 'ymax':1}
        def checkArr( A ) :
            if isinstance(A,(list,tuple)) : return A[0], A[1]
            else: return -A,A
                
        if 'xlim' in kargs : 
            domain['xmin'], domain['xmax'] = checkArr(kargs['xlim'])
        if 'ylim' in kargs : 
            domain['ymin'], domain['ymax'] = checkArr(kargs['ylim'])
                
        # -----------------------------------------
        if ( (xsamples is None and ysamples is None) and  
             (xplane   is None and yplane   is None) ) : return self

        numVerts = self.numSeg + 1
        vertSingle,indicesSingle = self._getSlice(g,numVerts,domain,xplane, yplane)
        vertsSet,indicesSet = self._getSliceSet(g,numVerts,domain,xsamples, ysamples)
        
        if vertSingle is not None :
            verts,indices = vertSingle,indicesSingle
            if vertsSet is not None :
                startIndex = len(verts)
                verts = np.concatenate( (verts,vertsSet), axis=0)
                shiftedIndicesSet =  [ [ i+startIndex for i in line ] for line in indicesSet] 
                numbIndiceSingle = len(indicesSingle)
                numbIndicesSet = len(shiftedIndicesSet)
                temp = [None]*(numbIndiceSingle + numbIndicesSet)
                for i in range(numbIndiceSingle) : temp[i] = indicesSingle[i]
                for i in range(numbIndicesSet) : temp[i+numbIndiceSingle] = shiftedIndicesSet[i]
                indices = temp

        else : 
            verts,indices = vertsSet,indicesSet

        self.lineIndices = indices
        self.vertexCoor = verts
        self.segmIndices = indices
        segs = self._get_segments( verts, indices  )
        self.set_segments(segs)
        if len(self.get_color() ) == 1 :
            c = self.get_color()[0]
            self.set_color( [c]*len(indices)   )

        self._set_geometric_bounds()
        self._planeNormal = planeNormal
        return self

    def map_xySlice_from_op(self, operation, **kargs ) :
        """
        Using a functional operation, lines of z=f(x) or z=f(y) in domain -1 to 1.

        Parameters
        ----------
        operation : function object 
            Function that takes one argument, a Numpy float
            3 X N array of xyz coordinates in the domain.
            The return value is a 3 X N array of xyz coordinates.

        xplane : contour at a constant value of x=xplane

        yplane : contour at a constant value of y=xplane

        xset :   number of evenly spaced contours, y-axis normal plane

        yset :   number of evenly spaced contours, x-axis normal plane

        xlim :   X domain of operation ( default: [-1,1] )

        ylim :   Y domain of operation ( default: [-1,1] )

        Returns
        -------
        self : ParametricLine object

        """

        def zValue_generator(xArr, yArr, op) :
            xlen = len(xArr)
            ylen = len(yArr)
            xcol = np.tile( [1.0,0.0,0.0], (xlen,1))
            xVal = np.multiply(xcol.T,xArr).T
            xVal3 = np.tile(xVal, (ylen,1))
            
            ycol = np.tile( [0.0,1.0,0.0] , (ylen,1))
            yVal = np.multiply(ycol.T,yArr).T
            yVal3 = np.repeat(yVal,xlen,axis=0)
            
            pts = xVal3 + yVal3
            pts = np.reshape(pts,(xlen*ylen,3))
            x,y,z = op(pts.T)
            z = np.reshape(z,(ylen,xlen))
            # the following, compensates for interp2d for datagrid
            # which drops the first dimension, see _getSlice 
            if ylen == 1 : z = z[0]
            return z

        def g(xArr, yArr) :
            return zValue_generator(xArr, yArr, operation)
        # -----------------------------------------
        name = _getFunctionName(operation,self._geomName)
        if name is not None : self._geomName = name
        return self._map_xySlice( g, **kargs  )

    def map_xySlice_from_datagrid(self, datagrid, scale=None, **kargs  ) :
        """
        Using a datagrid, lines of z=f(x) or z=f(y) in domain [-1, 1]

        Parameters
        ----------
        datagrid : 2D float array

        scale: multiplier for the Z direction.

        xplane : contour at a constant value of x=xplane

        yplane : contour at a constant value of y=xplane

        xset :   number of evenly spaced contours, y-axis normal plane

        yset :   number of evenly spaced contours, x-axis normal plane

        Returns
        -------
        self : ParametricLine object

        """
        if 'xlim' in kargs or 'ylim' in kargs :
            warnings.warn('xlim and ylim args for Datagrid xySlice not available.')
            return self
            
        data = datagrid
        dmax = np.amax(datagrid)
        dmin = np.amin(datagrid)
        delta = dmax-dmin
        data = (datagrid - dmin)/delta
        xd = np.linspace(-1, 1, data.shape[0] )
        yd = np.linspace(-1, 1, data.shape[1] )
        g =  interpolate.interp2d(yd, xd, data, kind='cubic')

        # -----------------------------------------
        self._map_xySlice( g, **kargs  )
        if scale is not None : 
            self.vertexCoor = np.multiply(self.vertexCoor,[1,1,scale])
            indices = self.segmIndices
            verts= self.vertexCoor
            segs = self._get_segments( verts, indices  )
            self.set_segments(segs)
            
        return self

    def scale_dataframe(self,X,Y,Z) : 
        """
        Scaling the line geometry based on a datagrid.

        Parameters
        ----------
        X, Y, Z : N x M arrays
            Minimum and maximum values of the arrays are used to scale
            and translate the surface from an intial domain of
            [ (-1,1), (-1,1), (0,1) ]

        Returns
        -------
        self : ParametricLine object

        """

        return _scale_dataframe(self,X,Y,Z)

class SegmentLine(ColorLine3DCollection) :
    
    def __init__( self, vertexCoor , name=None, **kwargs ) :
        """
        A single line of squential segments from an array of vertices.

        Parameters
        ----------
        vertexCoor : V x 3 numerical array-like. 
            An array of V number of xyz vertex coordinates.

        name : string, optional, default: None.
            Descriptive identifier.
                
        Other Parameters
        ----------------
        **kwargs
            All other parameters are passed on to 's3dlib.surface.ColorLine3DCollection'.

        """       
        line = vertexCoor  # <--- super will make this to a numpy float array
        indices = [ [*range(len(line))] ]
        super().__init__(line, indices, name, **kwargs)



# =================================================================================+
# DevNote : Although Surface, Vector and Line collection classes have similar      |
#           methods (eg. map..., name, clip ..), differences among these classes   |
#           are sufficiently different to not warrant a class that these three     |
#           classes could use with multiple inheritence.  Instead, these           |
#           methods are common among these three classes.                          |
#           FutDev may be useful to define such a class, but MROI.                 |
#                                                                                  |
#           Overall ColorLine3DCollection class needs a architectural change.      |
#           Possibly a list of ColorSegment3DCollection objects which are          |
#           contained within each ColorLine3DCollectio object. Each list           |
#           is a basically a line.  Analogous to Surface collections which         |
#           which are a collection of polyhedrons. This would be an advantage      |
#           for writing and reading to files.                                      |
#                                                                                  |
#           Becoming a priority: Need to have separate object properties for       |
#           face colors and edge colors, and don't use the Matplotlib properties   |
#           for these in Poly3DCollection, Line3DCollection. Too much 'magic'      |
#           happening which is out of control.  Then use set_color                 |
#           when needed.  This applies to both surface and line objects.           |
#           This also would allow constructing classes independent of              |
#           Poly3DCollection, Line3DCollection for future developement,            |
#           allowing an interface to matplotlib to be constructed instead          |
#           inheritance. So:                                                       |
#                                                                                  |
#           Still using Matplotlib 'private' properties, see:                      |
#               method:      _postProc_surfaceColors                               |
#               properties:  _dfacecolors and _dedgecolors                         |
#                                                                                  |
#           Unresolvable problems of using Matplotlib :                            |
#               1. can't independently set linewidth for each edge.                |
#                  Result: edges are shown for composites with different           |
#                  levels of transparency.                                         |
#               2. Poly3DCollection and Line3DCollection do not play               |
#                  nice on the same plot.  Result: the different class             |
#                  objects can't be added together resulting in inaccurate         |
#                  overlaying of faces and lines in a view.                        |
#                                                                                  |
# =================================================================================+

class KWprocessor():
    """  Access values from defaults or kwargs. """

    def __init__(self,valDict, methodStr=None, checkInputKeys=True, **kargs) :
        """
        Validate input of a method's kwargs.  Internal use ONLY.

        Parameters
        ----------
        valDict : dictionary of argument names, defaults and valid values.
            key -   keyword argument name
            value - single default value, or list of 2D lists to be validated:
                    If a list of list, the structure is:
                    value[0] - ['default value']    <- single value
                    value[1] - [ a, b, c, ..... ]   <- list of a list of valid values

        methodStr : string, optional, default: None
            Identifier of the method called where a validation error occured
            which is shown in the error message.

        checkInputKeys : bool, default: True
            If True, check if all kargs are in the valDict.
            Unknown keys will produce a ValueError.
            If False, ignore kargs not in valDict.

        kargs :
            The dictionary object containing the input keyword:values to be extracted.

        """
        # FutDev: could use a dictionary for value[1] to validate ranges and types.
        self.valDict = valDict
        self.kargs = kargs
        self.method = methodStr
        if checkInputKeys :
            if len(kargs) == 0    : return
            if len(valDict) == 0  : return
            errMsg = "Keyword, {}, is not recognized. \n List of possible keys are: {}"
            if methodStr is not None : errMsg = '['+methodStr+'] - ' + errMsg
            for key in kargs.keys() : 
                if not key in valDict.keys() :
                    raise ValueError( errMsg.format(key, list(valDict.keys())  ) )
        return    

    def _defVal(self,v) :
        value = self.valDict[v]
        if isinstance(value,list) :
            firstVal = value[0]
            return firstVal[0] if isinstance(firstVal,(list,tuple)) else value
        return value

    def _valueKeywordError(self,keyStr,usrVal) :
        dicVal = self.valDict[keyStr]
        errList = dicVal[1]
        for i in range(2,len(dicVal) ) : errList += dicVal[i] 
        errA = "Keyword, {}, value of {} is not recognized.\n List of possible values are: {}"
        errB = "Keyword, {}, value of '{}' is not recognized.\n List of possible values are: {}" 
        errStr =  errB if isinstance(usrVal,str) else errA      
        if self.method is not None : errStr = '['+self.method+'] - ' + errStr
        return errStr.format(keyStr,str(usrVal),errList)
        
    def getVal(self,keyStr,checkValidity=True)  :
        """
        Return the named argument value.

        Parameters
        ----------
        keyStr : string. Keyword

        checkValidity : bool, optional, default: True
            If True, karg value is checked to be valid
            among the list of valid values.

        Raises
        ------
        ValueError
            if invalid keyword assignment.

        Returns
        -------
            value of keyword, or the default value if the keyStr is not found in kargs.

        """
        def listOfList(v) :
            if not isinstance(v,(list,tuple)) : return False
            return isinstance(v[0],(list,tuple)) 
        kargs = self.kargs
        if len(kargs) == 0 : return self._defVal(keyStr)
        if keyStr not in self.valDict : # only for development purposes.
            print("################ CODE ERROR KW001 #####################")
            return None
        if keyStr in kargs : 
            dicVal = self.valDict[keyStr]
            usrVal = kargs[keyStr]
            if listOfList(dicVal) :
                for vLst in dicVal :
                    if usrVal in set(vLst) : return vLst[0] #<< value is the first.
                if not checkValidity : return usrVal
                else :
                    raise ValueError( self._valueKeywordError(keyStr,usrVal))
            return usrVal
        return self._defVal(keyStr)

    def filter(self):
        """
        Removes defaults from kargs.
        Use implies that more than default kargs are to be used,
        thus 'checkInputKeys' should be set to False at __init__.
        """
        for key in self.valDict:
            try :   del self.kargs[ key ]
            except: pass
        return self.kargs


# FutDev: include transform

def _axisUSF(axes) :
    # axis unit scale factor
    xlim = axes.get_xlim()
    ylim = axes.get_ylim()
    zlim = axes.get_zlim()
    xR = xlim[1]-xlim[0]
    yR = ylim[1]-ylim[0]
    zR = zlim[1]-zlim[0]
    sf = np.array([ yR*zR, xR*zR, xR*yR ])
    return sf/np.linalg.norm(sf)

def _coor_dotprod(direction,cIndex,refCoor) :
    # refCoor may be vertex(contours) or face(clip) coordinates.
    # determines the distance from the refCoor points to either:
    #    0. a plane, (planar coordinates)
    #    1. a line,  (cylindrical coordinates)
    #    2. a point, (spherical coordinates)
    # used for surface contours (vertices), surface clipping (faces)
    # and line clipping (vertices).
    verts = refCoor
    if cIndex == 1 :      # ..................... cylindrical
        if direction is None :
            dotprod = np.linalg.norm(verts[:,:2],axis=1)  # z-axis default
        else :
            unitDirection = np.divide( direction, np.linalg.norm(direction) )
            dmag = np.dot(verts,unitDirection)
            uVect = np.tile(unitDirection,(len(dmag),1) )
            dVect = np.array([dmag]).T
            zeta = verts - uVect*dVect
            dotprod = np.linalg.norm(zeta,axis=1) 
    elif cIndex == 2 :    # ..................... spherical
        dotprod = np.linalg.norm(verts,axis=1)           
    else :                # ..................... planar
        if direction is None : direction = [0,0,1.0]  # xy-plane default
        unitDirection = np.divide( direction, np.linalg.norm(direction) )
        dotprod = np.dot(verts,unitDirection)
    return dotprod

def _getFunctionName(op,xName=None) :
    if xName is not None : return xName
    # assign function name if not a lambda function
    fobar = lambda x : x
    if op.__name__ != fobar.__name__ :
        name = op.__name__
    else :
        name = None
    return name

def _set_geomBounds(v,bounds) :
    """Set _bounds dictionary values from vertex coordinates."""

    # v = self.vertexCoor
    xlim,ylim,zlim = np.transpose( [np.amin(v ,axis=0), np.amax(v ,axis=0) ] )
    bounds['xlim'] = xlim
    bounds['ylim'] = ylim
    bounds['zlim'] = zlim
    rmax_xy =  np.amax( np.linalg.norm(v[:,0:2], axis=1) )
    rmin_xy =  np.amin( np.linalg.norm(v[:,0:2], axis=1) )
    bounds['r_xy'] = [rmin_xy, rmax_xy ] 
    rmax_xyz = np.amax( np.linalg.norm(v, axis=1) )
    rmin_xyz = np.amin( np.linalg.norm(v, axis=1) )
    bounds['rorg'] = [rmin_xyz, rmax_xyz ] 
    bounds['rlim'] = [rmax_xy, rmax_xyz ]  # .. vestigal
    return    

def _scale_dataframe(obj,X,Y,Z) :
    """
    Linear scale and translate the geometry.

    Used for scaling a geometry based on a datagrid.

    Parameters
    ----------
    X, Y, Z : 2d arrays
        Minimum and maximum values of the arrays are used to scale
        and translate the surface from an intial domain of
        [ (-1,1), (-1,1), (0,1) ]

    Returns
    -------
    obj :  first calling argument, object 

    """

    Xmin, Xmax = np.amin(X), np.amax(X)
    Ymin, Ymax = np.amin(Y), np.amax(Y)
    Zmin, Zmax = np.amin(Z), np.amax(Z)
    Xscale = (Xmax-Xmin)/2
    X0 = (Xmax+Xmin)/2
    Yscale = (Ymax-Ymin)/2
    Y0 = (Ymax+Ymin)/2
    Zscale = Zmax-Zmin
    Z0 = Zmin
    obj.transform(scale=[Xscale,Yscale,Zscale],translate=[X0,Y0,Z0])
    return obj

def _transformAxis(selfProp, lenmult=1.0, width=1.5, colors=None, negaxis=True) :
    """
    Line3DCollection of a transformed coordinate axis.
    
    just take method from surface version 1.0 and allow use for all
    other transformable objects.
    Called from the interface method get_transformAxis
    
    Returns
    -------
    Line3DCollection

    """
    self_scale = selfProp['scale']
    self_rotation = selfProp['rotation']
    self_translation = selfProp['translation']
    scale = np.multiply(lenmult,self_scale)
    dDelta = np.multiply(np.identity(3).tolist(),scale).tolist()
    npdot = np.dot(dDelta,self_rotation)
    npdot2 = np.expand_dims(npdot,axis=1)
    if negaxis :
        lines = np.concatenate( (npdot2,-npdot2),axis=1)
    else :
        lines = np.insert(npdot2,0,[0,0,0],axis=1)
    lines = np.add(lines,self_translation)
    if colors is None : 
        colorMap = ['r','g','b']
    else :
        colorMap = colors
    axisCol = Line3DCollection(lines, colors=colorMap, linewidths=width)
    return axisCol      

def _normLength( N, exp, fac, limitarea=None) :
    """
    Default face normal vector length.
    
    Length proportional to the mean triangulated face edge length.
    """
    if limitarea is None : limitarea = 4*np.pi
    # aratio - area ratio as a function of the number of faces, N.
    # heuristic values, exp and fac, are adjustment for small rez surfaces,
    # with aratio limit to 1 for large N. (continuous smooth surface)
    aratio = 1.0 - fac*np.power(N,-exp)
    sc = 0.5*(1.520)  # scale of normalized edge length
    nLen = sc*np.sqrt(aratio*limitarea/N)
    return nLen

def _interpret_domain(domain) :
    """
    Process domain argument for implsurf and CubicSurface.domain to
    the domain for x,y,z domains.
    """
    # _interpret_domain no longer a static method... update 1.3.0 from 1.2.0
    errCode = 0
    if isinstance(domain, (int,float)) : 
        dftDm = np.array( [ [-domain,domain ], [-domain,domain ], [-domain,domain ] ] )
    else :
        correctDomain = True
        try :
                npDm = np.array(domain, dtype=float)
                shDm = npDm.shape
                nbDm = npDm.ndim
                correctDomain = nbDm==1 or nbDm==2
                errCode = 1
        except :
            correctDomain = False
            errCode = 2
        if correctDomain :
            if nbDm == 1 : 
                dftDm = np.array( [ domain, domain, domain] )
                correctDomain = shDm[0] == 2
                errCode=3
            elif nbDm == 2 : 
                dftDm = np.array( domain )
                correctDomain = shDm[0]==3 and shDm[1]==2
                errCode=4
            else :
                correctDomain = False
                errCode=5
        if not correctDomain :
            raise ValueError("Incorrect domain argument passed to impf.", errCode)
    return dftDm


# =================================================================================+
# Utility Functions:                                                               |
# Useful for construction of 3D Matplotlib figures composed of S3Dlib objects.     |
# =================================================================================+

def setupAxis(axes, **kargs ) :
    """
    Add XYZ origin coordinate axis Vector3DCollection to the axes.

    Parameters
    ----------
    axes : Matplotlib 3D axes object.

    length  : coordinate axis length (single or 3-value list), default: 1.5

    color   : axis color (single or 3-value list), default: 'black'

    offset  : axis vector offset from the origin (single or 3-value list), default: 0

    labels  : axis label (single or 3-value list), default: ['X','Y','Z']

    pad     : label padding from the vector axis head (single or 3-value list), default: 0.15

    width   : axis line widths, default: 2

    negaxis : show negative axis line, default: False

    alr     : axis length ratio, head size to axis length, default: 0.2

    Returns
    -------
    axes

    """

    kargDft = { 'length': 1.5, 'color':'k', 'labels':['X','Y','Z'],
                'offset':0.0, 'pad': 0.15, 'width':2, 'alr':0.2, 'negaxis':False}
    KW = KWprocessor(kargDft,'setupAxis',**kargs)
    length = KW.getVal('length')
    color = KW.getVal('color')
    offset = KW.getVal('offset')
    labels = KW.getVal('labels')
    pad = KW.getVal('pad')
    width = KW.getVal('width')
    negaxis = KW.getVal('negaxis')
    alr = KW.getVal('alr')
    
    # Note: not checking array length==3, or a valid color.
    if color is None : color='k'
    if not colors.is_color_like(color) :  # check if arg is a color array, not an array of colors 
        color = np.array ( [colors.to_rgba(color[0]), colors.to_rgba(color[1]), colors.to_rgba(color[2])  ] )
    else :
        color = np.array ( [colors.to_rgba(color), colors.to_rgba(color), colors.to_rgba(color)]  )

    if not isinstance(length,(list,tuple,np.ndarray)) : length =  [length,length,length] 
    length=np.array(length)

    if not isinstance(offset,(list,tuple,np.ndarray)) :  offset = [offset,offset,offset]
    offset = np.array(offset)

    if labels is None : labels = ['X','Y','Z']
    if not isinstance(labels,(list,tuple,np.ndarray)) :  labels = [labels,labels,labels]
    
    if not isinstance(pad,(list,tuple,np.ndarray)) :  pad = [pad,pad,pad]
    pad = np.array(pad)

    """Set XYZ origin coordinate axis arrows."""
    seglength = length-offset
    endCoor = np.multiply(np.identity(3),seglength)
    strCoor = np.multiply(np.identity(3),offset)

    xyz = np.transpose( strCoor )
    uvw = np.transpose( endCoor )
    vcf = Vector3DCollection(xyz,uvw,alr=alr, color=color[0], linewidth=width)
    vcf._set_vectColor(color)
    axes.add_collection3d( vcf )

    if negaxis :
        endCoor = np.multiply(np.identity(3), -1.0*length)
        strCoor = np.multiply(np.identity(3), -1.0*offset)
        xyz = np.transpose( strCoor )
        uvw = np.transpose( endCoor )
        vertexCoor = np.concatenate( (endCoor,strCoor),axis=0 )
        segmIndices = [ [0,3],[1,4],[2,5] ]
        line = ColorLine3DCollection(vertexCoor,segmIndices, color=color, linewidth=width)
        axes.add_collection3d( line )

    pos = length + pad
    axes.text(pos[0],0,0,labels[0], color = color[0], fontweight='bold', fontsize='large',
            horizontalalignment='center', verticalalignment='center')
    axes.text(0,pos[1],0,labels[1], color = color[1], fontweight='bold', fontsize='large',
            horizontalalignment='center', verticalalignment='center')
    axes.text(0,0,pos[2],labels[2], color = color[2], fontweight='bold', fontsize='large',
            horizontalalignment='center', verticalalignment='center')
    return axes     

def standardAxis( axes, **kargs ) :
    """
    Add XYZ origin coordinate axis Vector3DCollection and
    set the axis viewing angle, projection type and clear planes.
    
    Parameters
    ----------
    axes : Matplotlib 3D axes object.

    kargs: identical to those for setupAxis function.

    Returns
    -------
    axes
    
    """
    setupAxis(axes, **kargs)
    minmax = (-1,1)
    axes.set(xlim=minmax, ylim=minmax, zlim=minmax)
    axes.view_init(30, 30)
    axes.set_proj_type('ortho')
    axes.set_axis_off()
    return axes

def auto_scale(axes,*obj3d,**kargs) :
    """
    Scale Axes3d using object bounds.

    Parameters
    ----------
    axes : Matplotlib 3D axes.

    obj3d : positional arguments of s3dlib objects
    
    uscale : scaling in the range (0.5,2.5)
        if neither uscale or rscale are defined - 
        .    Each axes ranges are independently
        .    set by sizes of input objects.
        if both uscale and rscale are defined -
        .    uscale will be ignored.
        if defined using an int or float -
        .    Each axes ranges have the same scalar
        .    length, set by the values of the objects.
        .    Out of range values will set the scaling
        .    to the minimum or maximum.
        else - scaling is one.
    
    rscale : scaling in the range (0.5,2.5)
        if neither uscale or rscale are defined - 
        .    Each axes ranges are independently
        .    set by sizes of input objects.
        if both uscale and rscale are defined -
        .    uscale will be ignored.
        if defined using an int or float -
        .    Each axes ranges have the same minimum
        .    and maximum, set by the values of the
        .    objects radial extent from the origin.
        .    Out of range values will set the scaling
        .    to the maximum.
        else - scaling is one.

    Returns
    -------
    axes, which has been scaled

    """   
    X,Y,Z,R = [],[],[],[]
    for obj in obj3d :
        X.append(obj.bounds['xlim'][0])
        X.append(obj.bounds['xlim'][1])
        Y.append(obj.bounds['ylim'][0])
        Y.append(obj.bounds['ylim'][1])
        Z.append(obj.bounds['zlim'][0])
        Z.append(obj.bounds['zlim'][1])
        R.append(obj.bounds['rorg'][1])
    
    if 'uscale' in kargs :
        sc = kargs['uscale']
        if not isinstance(sc,(int,float)): sc=1.0
        if sc > 2.0 : sc = 2.0
        if sc < 0.5 : sc = 0.5
        AxLen = lambda a : max(a)-min(a)
        AxMpt = lambda a : 0.5*( max(a) + min(a) )
        maxAlen = sc*max( [AxLen(X), AxLen(Y), AxLen(Z)]  ) 
        X = [ (AxMpt(X)-0.5*maxAlen), (AxMpt(X)+0.5*maxAlen) ]
        Y = [ (AxMpt(Y)-0.5*maxAlen), (AxMpt(Y)+0.5*maxAlen) ]
        Z = [ (AxMpt(Z)-0.5*maxAlen), (AxMpt(Z)+0.5*maxAlen) ]

    if 'rscale' in kargs :
        sc = kargs['rscale']
        if not isinstance(sc,(int,float)): sc=1.0
        if sc > 2.0 : sc = 2.0
        if sc < 0.5 : sc = 0.5
        rmax = sc*max(R)
        minmax = (-rmax,rmax)
        axes.set(xlim=minmax, ylim=minmax, zlim=minmax )
        return axes      

    axes.auto_scale_xyz(X,Y,Z)
    return axes

def eulerRot(theta, phi, psi=0, useXconv=True, inrad=False) :
    """Rotation matrix from Euler angles, X or Y convention."""
    # X-convention (zxz') and Y-convention (zyz') rotations.
    tht_Rad = theta
    phi_Rad =  phi
    psi_Rad =  psi
    if not inrad :
        tht_Rad = theta*np.pi /180.0
        phi_Rad =  phi*np.pi /180.0
        psi_Rad =  psi*np.pi /180.0

    cT,sT = np.cos(tht_Rad),  np.sin(tht_Rad)
    cP,sP = np.cos(phi_Rad),  np.sin(phi_Rad)
    cS,sS = np.cos(psi_Rad),  np.sin(psi_Rad)
    T_t = [ [ cT, sT, 0 ], [-sT, cT, 0 ], [ 0 , 0 , 1 ]  ]
    T_X = [ [ 1 , 0 , 0 ], [ 0 , cP, sP], [ 0 ,-sP, cP]  ]
    T_Y = [ [ cP, 0 ,-sP], [ 0 , 1 , 0 ], [ sP, 0 , cP]  ]
    T_s = [ [ cS, sS, 0 ], [-sS, cS, 0 ], [ 0 , 0 , 1 ]  ]
    if useXconv : T_p = T_X
    else:         T_p = T_Y
    M = np.matmul( T_p, T_t )
    if psi != 0 : M = np.matmul( T_s, M )
    return M

def axisRot(alpha, direction, inrad=False) :
    """Rotation matrix from rotation about a directional axis"""
    # see Wikipedia: Quaternions and spatial rotation
    alp_Rad = -alpha
    if not inrad :
        alp_Rad = -alpha*np.pi /180.0

    u = np.divide( direction, np.linalg.norm(direction) )
    sT = np.sin(alp_Rad/2.0)
    qR = np.cos(alp_Rad/2.0)
    qI, qJ, qK = sT*u[0], sT*u[1], sT*u[2] 
    M = [
        [ (1-2*(qJ*qJ + qK*qK)),    2*(qI*qJ - qK*qR) ,    2*(qI*qK + qJ*qR)  ] ,
        [    2*(qI*qJ + qK*qR) , (1-2*(qI*qI + qK*qK)),    2*(qJ*qK - qI*qR)  ] ,
        [    2*(qI*qK - qJ*qR) ,    2*(qJ*qK + qI*qR) , (1-2*(qI*qI + qJ*qJ)) ]
    ]
    return M

def vectRot(xDirection, in_planeDirection) :
    """Rotations matrix from two vectors """
    ipDir = np.array(in_planeDirection)
    xDir = np.array(xDirection)
    zDir = np.cross(xDir,ipDir)
    yDir = np.cross(zDir,xDir)
    vect = np.array( [xDir,yDir,zDir] )
    vnorm = np.linalg.norm(vect,axis=1)
    unitVects = np.divide(vect.T,vnorm).T
    return unitVects

def rtv(direction,elev=None,azim=None) :
    """Transform direction vector 'relative to view' axis viewing angles."""
    if elev is None : elev = _DFT_VIEW[0]
    if azim is None : azim = _DFT_VIEW[1]
    azim_Rad = azim*np.pi /180.0
    elev_Rad = elev*np.pi /180.0
    cA,sA = np.cos(azim_Rad),  np.sin(azim_Rad)
    cE,sE = np.cos(elev_Rad),  np.sin(elev_Rad)
    TM = [ [ cE*cA, cE*sA, sE ],  [ -sA, cA, 0.0 ],  [  -sE*cA, -sE*sA, cE ] ]
    reldir = np.dot(direction,TM)
    return reldir

def elev_azim_2vector(elev, azim, inrad=False) :
    """Unit vector in the terms of angular position."""
    tht_Rad = azim
    phi_Rad = np.pi/2 - elev
    if not inrad :
        tht_Rad = azim*np.pi /180.0
        phi_Rad = (90.0-elev)*np.pi /180.0
    k = np.cos(phi_Rad)
    r = np.sin(phi_Rad)
    i = r*np.cos(tht_Rad)
    j = r*np.sin(tht_Rad)
    return np.array( [i,j,k] )

def frame_to_value(x, *valList, **kargs) :
    """
    Y = f(x) from a collection of linear segment functions f(x)
    
    Parameters
    ----------
    x : float value between 0 and 1.

    valList : list of values, evenly spaced along x. For multiple lists, all must be the same length.

    stepped : boolean, default is False.
        If True, function evaluate to a constant of the first value between
        valList values.  If False, function is linear between valList values.
        For True, number of linear segments is the valList length.
        For False, number of linear segments is one less than valList length.

    Returns
    -------
    List of values of size valList arguments number.

    """
    # 0 < x < 1
    def vFromList(f,A) :
        nSeg = len(A) - 1
        i = math.floor(f*nSeg)
        if i>=nSeg : return A[nSeg]
        x = (f*nSeg) % 1
        val = A[i] + x*( A[i+1]-A[i] )
        return val
    def vFromStep(f,A) :
        nSeg = len(A)
        i = math.floor(f*nSeg)
        if i>=nSeg : return A[nSeg-1]
        if i<0 : return A[0]
        val = A[i]
        return val
    stepped = False
    if 'step' in kargs : stepped = kargs['step']   
    vals = []
    for vList in valList :
        if stepped : vals.append( vFromStep(x,vList) )
        else       : vals.append( vFromList(x,vList) )
    return vals

def density_function(vals, bins=(10,10), scale=True, xbnd=True, ybnd=True, kind='linear') :
    """
    2D histogram density, Data_density = f(x,y) from a x,y dataset.
    
    Parameters
    ----------
    vals : 2 X N array of x,y data values.

    bins : [int,int], the number of bins for each dimension.

    scale : bool or number, default: True
        If True, normalize values to the average.
        If False, scale is 1, the density value of data per
        bin area.  When numerical, multiplier for
        the density value.

    xbnd : bool, default : True
        If True, x values are extended for density of zero.
        For cyclic mapping, this should be set to False.

    ybnd : bool, default : True
        If True, y values are extended for density of zero.
        For cyclic mapping, this should be set to False.

    kind : {‘linear’, ‘cubic’, ‘quintic’}, default: 'linear'
        The kind of spline interpolation to use.
   
    Returns
    -------
    A function f(x,y)

    """
    vals = np.array(vals)
    bins = np.array(bins)
    Zvals,xedges,yedges = np.histogram2d(*vals,bins=bins)
    X,Y = 0,1  
    minval = np.amin(vals,axis=1)
    maxval = np.amax(vals,axis=1)
    datarng = (maxval-minval)
    area = datarng[0]*datarng[1]
    delta = datarng/bins
    endval = maxval + delta/2
    if xbnd :
        zeros = np.zeros( (Zvals.shape[X]+2,Zvals.shape[Y]) )
        zeros[1:bins[X]+1,:] = Zvals
        Zvals = zeros
        xedges = np.array([*(xedges-delta[X]/2),endval[X]])
    else :
        xedges = xedges+delta[X]/2
        xedges = xedges[:len(xedges)-1]
    if ybnd :
        zeros = np.zeros( (Zvals.shape[X],Zvals.shape[Y]+2) )
        zeros[:,1:bins[Y]+1] = Zvals
        Zvals = zeros
        yedges = np.array([*(yedges-delta[Y]/2),endval[Y]])
    else :
        yedges = yedges+delta[Y]/2
        yedges = yedges[:len(yedges)-1]
    
    if isinstance(scale,bool) : 
        scale = area*(bins[0]*bins[1])/len(vals[0]) if scale else 1.0

    Zvals = scale*Zvals/area
    f = interpolate.interp2d(yedges,xedges,Zvals,kind=kind)

    def density(xa,ya):
        size = min( len(xa),len(ya) )
        vals = [ f(ya[i],xa[i])[0] for i in range(size) ]
        return np.array(vals)

    return density

def get_points_from_cloud(cloud,domain=1,cmap=None,fade=None) :
    """
    xyz coordinate points, marker colors array, and ScalarMappable object.
    
    Parameters
    ----------
    cloud : a (M,N,P) shaped array of point cloud values.

    domain : a number, list or array, default: 1
        The domain of the function evaluation.  For a number, n,
        the x,y,z axes domains will be [-n,n]. For a 1-dimensional 2-element list,
        [a,b] will be assigned the domain for all 3 axes.  Using a list of list (array),
        as [ [a,b],[c,d],[e,f] ], the domain is assigned individually for each of the
        three coordinate axes.

    cmap : str or Colormap, default: Matplotlib default colormap
        A Colormap instance or registered colormap name.
        
    fade : an float, default: None
        Powerlaw adjust point opacity for emphasizing upper/lower values.
        The absolute values of fade must be in the range 0.1 to 10.

    Returns
    -------
    xyz : a 3 X S array, where S is the number of points in the domain.

    colors : array of point colors.

    ScMap : a Matplotlib ScalarMappable object.

    """
    dftDm = _interpret_domain(domain)
    objMap = _DFT_CMAP if cmap is None else cmap
    if isinstance(objMap,str) : 
        objMap = mpl.colormaps[objMap]
    S = np.array(cloud.shape)
    rngX,rngY,rngZ = dftDm[0],dftDm[1],dftDm[2]
    xyz = np.mgrid[ rngX[0]:rngX[1]:S[0]*1j,
                    rngY[0]:rngY[1]:S[1]*1j,
                    rngZ[0]:rngZ[1]:S[2]*1j ]
    sxyz = np.reshape(xyz,(3,-1))
    vxyz = cloud.flatten()
    normalizedData = (vxyz-np.min(vxyz))/(np.max(vxyz)-np.min(vxyz))
    vcolors = objMap(normalizedData)
    if fade is not None:
        if abs(fade) < 0.1 or abs(fade) > 10.0 :
            raise ValueError('[get_points_from_cloud]: fade out of range, 0.1 >= abs(fade) <= 10 ' )
        if fade > 0 :
            vcolors[:,3] = normalizedData**fade    # set color with gradient transparency.
        else :
            vcolors[:,3] = (1-normalizedData)**(-fade)
    norm = colors.Normalize(vmin=np.min(vxyz),vmax=np.max(vxyz))
    sm = cm.ScalarMappable(cmap=objMap, norm=norm)
    sm.set_array([])
    return sxyz, vcolors, sm


# =================================================================================+
# Only vertices and face indices are exported/imported.                            |  
# FutDev, include :                                                                |
#   1. metadata                                                                    |
#   2. face/vertex color                                                           |
#   3. vertex normals for external use for Gouraud shading, etc.                   |
#   4. line objects                                                                |
#   5. code in separate eximport module                                            |
#   6. file types (ascii,binary) : stl, ply, X3D, etc.                             |
# =================================================================================+

def get_surfgeom_from_obj(file, **kargs) :
    """ Create surface from obj formatted 'file'.
        kargs passed to instantiate surface in constructor. 
        'Proof of concept' for future versions.
    """
    # .....................................................................
    def file_records(in_file):
        try :
            f = open(in_file)
            buf = f.read()
            f.close()
        except :
            raise ValueError('ERROR: OBJ file not accessed !!!')
        for b in buf.split('\n'):
            b = ' '.join(b.split()) # remove extra spaces.
            if b.startswith('v '):
                yield ['v', [float(x) for x in b.split(" ")[1:]]]
            elif b.startswith('f '):
                datalist = [x for x in b.split(" ")[1:] ]
                # only get face indices at the first index...
                yield [ 'f', [  int(y.split("/")[0])-1   for y in datalist]   ]
            else:
                yield ['', ""]   
    def get_v_f(in_file):
        vertices,faces = [],[]
        for k, v in file_records(in_file):
            if   k == 'v': vertices.append(v)
            elif k == 'f': faces.append(v)
        if not len(faces) or not len(vertices):  return None, None
        return vertices, faces
    # .....................................................................
    v,f = get_v_f(file)
    v = np.array(v)
    v = v[:,[0,2,1]]
    v = np.multiply(v,[1,-1,1])
    if v is not None:  surface = Surface3DCollection(v,f,**kargs)
    else:              surface = SphericalSurface(0,'tetra',color='r')
    return surface

def save_surfgeom_to_obj(file,surface) :
    """ Export 'surface' geometry to obj formatted file. 
        'Proof of concept' for future versions.
    """
    v = surface.vertexCoor
    v = v[:,[0,2,1]]
    v = np.multiply(v,[1,1,-1])
    f = surface.fvIndices+1  # obj file indices start at 1
    vs = [ 'v '+np.array2string(x,prefix='',suffix='')[1:-1] for x in v ]
    fs = [ 'f '+np.array2string(x,prefix='',suffix='')[1:-1] for x in f ]
    geomLines = vs + fs
    with open(file,"w" ) as fileout :
        for line in geomLines :
            fileout.write(line+'\n')
        fileout.close()
    return
