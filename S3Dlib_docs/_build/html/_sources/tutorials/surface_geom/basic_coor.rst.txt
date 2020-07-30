.. _basic_coor_tut:

.. role::  raw-html(raw)
    :format: html



*****************************
Coordinates and Functions
*****************************

There are four basic surface classes, each with a native coordinate system.
These surfaces provide a basic grid of points which can be mapped using functional, image
or data methods.

Mapping for any surface is done in the respective native coordinates:

=======================    ==============================        ===================
Class                      Native Coordinates                    Geometric Shape
=======================    ==============================        ===================
PlanarSurface              x,y,z                                 square plate
PolarSurface               r, :raw-html:`&theta;`, z             circular disk
CylindricalSurface         r, :raw-html:`&theta;`, z             cylinder
SphericalSurface           r, :raw-html:`&theta; , &phi;`        sphere
=======================    ==============================        ===================

Surfaces are normalized to be contained within an x,y,z cube, with
center at the (0,0,0) and extending from -1 to 1 on each axis.  Scaling can be performed
either through the mapping methods or by a transform class method.

.. note::

   The four surface classes are derived from the same base class and use the same
   mapping method names with the same arguments, only differing in the native coordinates
   which are used in function definitions. 

.. note::

   Arguments passed to user-defined functions and their return values are Numpy arrays.




Function Coordinates
=========================================================================================

In the previous tutorial, the example surface was of type PlanarSurface.  As a result,
the function used to define the surface geometry used x,y,z Cartesian coordinates.
The return value of that function was a point in the same coordinates.
This function was::

    def planarfunc(xyz) :
        x,y,z = xyz
        r = np.sqrt( x**2 + y**2)
        Z = np.sin( 6.0*r )/2
        return x,y,Z

To use polar coordinates, a PolarSurface object is used with a function defined in
polar coordinates as:

.. literalinclude:: source/hw_polar.py
   :language: python
   :linenos:
   :emphasize-lines: 14,15

The only changes from the previous example is the use of a PolarSurface, line 14, and to use
the *polarfunc* function for the mapping operation, line 15.  The resulting plot is:

.. image:: images/hw_polar.png


.. _planar_coor_tut:

Planar
-----------------------------------------------------------------------------------------

The coordinate system for the PlanarSurface object is shown below. The
initial surface coordinates before mapping all have z=0.

.. image:: images/tut_plate.png
   :class: sphx-glr-single-img

For a PlanarSurface, the geometric mapping function will have the form:

.. literalinclude:: source/tut_codefrag_1.py
   :language: python
   :lines: 1-6




.. _polar_coor_tut:

Polar
-----------------------------------------------------------------------------------------

The coordinate system for the PolarSurface object is shown below. The
initial surface coordinates before mapping all have z=0 with
the radial coordinate domain :raw-html:`0 &le; r &le; 1`.



.. image:: images/tut_disk.png
   :class: sphx-glr-single-img

For a PolarSurface, the geometric mapping function will have the form:


.. literalinclude:: source/tut_codefrag_1.py
   :language: python
   :lines: 11-16



.. _cylindrical_coor_tut:

Cylindrical
-----------------------------------------------------------------------------------------

The coordinate system for the CylindricalSurface object is shown below. The
initial surface coordinates before mapping all have r = 1 with
the angular coordinate domain :raw-html:`0 &le; &theta; &lt; 2&pi;` and the vertical
domain :raw-html:`0 &le; z &le; 1`.



.. image:: images/tut_cylinder.png
   :class: sphx-glr-single-img

For a CylindricalSurface, the geometric mapping function will have the form:


.. literalinclude:: source/tut_codefrag_1.py
   :language: python
   :lines: 21-26




.. _spherical_coor_tut:

Spherical
-----------------------------------------------------------------------------------------

The coordinate system for the SphericalSurface object is shown below. The
initial surface coordinates before mapping all have r = 1 with
the angular coordinate domain :raw-html:`0 &le; &theta; &lt; 2&pi;` and the vertical
domain :raw-html:`0 &le; &phi; &le; &pi;`.

.. image:: images/tut_sphere.png
   :class: sphx-glr-single-img

For a SphericalSurface, the geometric mapping function will have the form:


.. literalinclude:: source/tut_codefrag_1.py
   :language: python
   :lines: 31-36







Parametric Functions
=========================================================================================


To exemplify plotting parametric functions, first consider a 2D line plot of
a simple power function based on the
`Simple Plot <https://matplotlib.org/tutorials/introductory/usage.html#the-object-oriented-interface-and-the-pyplot-interface>`_
example from matplotlib.


.. image:: images/tut_simple_plot.png
   :class: sphx-glr-single-img


Write the script in a three step process: define the function, create the 'pseudo' line objects, 'add' the lines
to the plot axis.


.. literalinclude:: source/tut_simple_plot.py
   :language: python


For 3D surfaces, use an analogous function definition:  the first argument is a coordinate location, the second
argument is the function parameter.

.. literalinclude:: source/tut_codefrag_1.py
   :language: python
   :lines: 41-47


So, to plot a parametric function of the previous 3D example using the three step process:

.. literalinclude:: source/tut_parametric_3D.py
   :language: python
   :emphasize-lines: 7,16,18,20

Which produces:

.. image:: images/tut_parametric_3D.png
   :class: sphx-glr-single-img

This example demonstrates the object-oriented approach to 3D surface plotting.  In Matplotlib, the
axis is an object.  With S3Dlib, surfaces are also considered objects and can be operated on
separately from the axis object.  In this example, one surface object was created from three 
surface objects through simple addition.  Then, that single surface object was added to the axis object.


This method of using a parametric function is demonstrated in numerous :ref:`example`.
The :ref:`dipole` example uses multiple surfaces on a single 3D axis.  The :ref:`param_set` example
plots each surface on an individual 3D axis.
