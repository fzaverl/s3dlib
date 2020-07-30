.. _intermediate:


****************************
Intermediate Surfaces
****************************

.. image:: images/intermediate.png
   :class: sphx-glr-single-img


It may be the case for complex surfaces to first construct an intermediate surface through mapping,
then apply additional functional mapping.  This technique was used in numerous examples.
This set of surfaces demonstrates the construction of simple intermediate surfaces.

The first set of surfaces on the top row were all constructed by geometric mapping a
*CylindricalSurface* object.  As a result, the surfaces all have cylindrical
coordinates as their native coordinate system.

The second set of surfaces on the bottom row were all constructed using the base class
*Surface3DCollection*. The object was instantiated by passing vertex coordinates and face indices
to the base class constructor.  This was followed by increasing the surface rez using
the *triangulate* method.  The coordinates and indices can either be directly constructed,
as was the case for the 'cube', or accessed from the subclass SphericalSurface,
as was the case for the 'icosahedron' and 'dodecahedron'.
Since these surfaces are of class *Surface3DCollection*, the native coordinates are xyz coordinates.

Note that when a colormap is applied to all surfaces in the highlighted line, the
z coordinate was used since this is the same for all the native coordinates for
these different objects.

.. literalinclude:: source/ex_intermediate.py
   :language: python
   :emphasize-lines: 114


