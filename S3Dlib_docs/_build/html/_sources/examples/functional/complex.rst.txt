.. _complex:

.. role::  raw-html(raw)
    :format: html


***********************************************************
Complex Number Representation, Geometry and Colormap
***********************************************************

Polar Coordinates
========================================================

.. image:: images/complex.png
   :class: sphx-glr-single-img


In this example, both functions are expressed in the native polar coordinates, z = 
:raw-html:`f(r,&theta;)`.
For the square root function, the real component is geometrically visualized 
whereas the imaginary component is expressed as color.  The square function in
the right plot expresses both real and imaginary components geometrically with
two surfaces combined as one.  Each function uses the *isReal* argument to
indicate whether the real or imaginary component is returned.

In polar coordinates, the domain for the angular coordinate is
:raw-html:`0 &le; &theta; &le; 2&pi;`.
The domain for the square root function is
:raw-html:`0 &le; &theta; &le; 4&pi;`,
hence the need for the highlighted statement in the following
code.  



.. literalinclude:: source/ex_complex.py
   :language: python
   :emphasize-lines: 18


.. _complex_planar:

Planar Coordinates
========================================================

.. image:: images/complex_2.png
   :class: sphx-glr-single-img

In this example, both functions are expressed in the planar coordinates, z = f(x,y).
As a result, the function definitions for this case are in x,y coordinates.

From de Moivre's formula, there is a positive and negative solution for
the square root.  Each solution requires a separate surface and therefore
two surfaces, positive and negative, were combined for the square root function.
The *isPos* argument sets whether the positive or negative surface 
coordinates are returned from the function.

.. literalinclude:: source/ex_complex_2.py
   :language: python
   :lines: 11-54
   :emphasize-lines: 6-13,15-19,25,29,33,37,39



Polar returning Planar
========================================================

For functions in one coordinate system, but return coordinates
in planar coordinates, the geometric mapping can return xyz
coordinates using the *returnxyz* argument instead of native
coordinates.  When this is the case, the function must
evaluate all x, y, and z coordinates.  For this example

| :raw-html:`x = g(r,&theta;)`
| :raw-html:`y = h(r,&theta;)`
| :raw-html:`z = f(r,&theta;)`

Although not needed in the current example, this
can easily be done by slightly modifying the functional
operations, as shown in the highlighted code below.
The  :ref:`dini` for the Dini surface is a good
example where returning x,y,z coordinates is applicable.

Using the following code, the resulting plots will be identical to those shown in
the first example.

.. literalinclude:: source/ex_complex_xyz.py
   :language: python
   :lines: 11-49
   :emphasize-lines: 11-13,19-21,27,29,32-35



Cylindrical Coordinates
========================================================


.. image:: images/complex_3.png
   :class: sphx-glr-single-img

In cylindrical coordinates,
( :raw-html:`r,&theta;,z` )
, there are two independent variables, 
:raw-html:`&theta; and z`, for the CylindricalSurface object.
The predefined surface object uses a constant radial coordinate for
this geometry with r = 1 and
:raw-html:`-1 &le; Z &le; 1`.
To apply the cylindrical surface to polar coordinates, use a simple 
linear tranformation for :raw-html:`0 &le; r &le; 1`, as

| :raw-html:`r =  &half; ( 1 - z )`.

where now the function for z can be expressed in terms of r 
and :raw-html:`&theta;`.
The r(z) linear function has a negative slope so that the surface
normals have a positve upward direction on the polar surface.

Compared to the polar coordinate plot, only a minor change was
needed in the function definitions. The only
noticeable difference in polar versus the cylindrical surface
is the minor resolution change in the square root plot.

.. note::
   This example is for demonstration purposes and cylindrical coordinates
   are not convenient to use for such a simple case.  However, it does
   demonstrate using a base surface coordinate with a range (in this case z)
   as a substitute for a coordinate which is constant (in this case r).
   When applied in this manner, coordinate transformations are necessary.

   The selection of using various base objects by considering the distribution
   of vertices in the native coordinates is further discussed in
   the :ref:`base_selection` guide.
   

.. literalinclude:: source/ex_complex_3.py
   :language: python
   :lines: 11-47
   :emphasize-lines: 8,16,25,30


Alternative Representations
========================================================

The two complex functions in this example have used two
different representations.  That being combining real and imaginary components
into one color mapped surface, or representing the two components
as two surfaces.  These representation may be used for either function,
as shown below by switching the representations.


.. image:: images/complex_4.png
   :class: sphx-glr-single-img


.. literalinclude:: source/ex_complex_4.py
   :language: python
   :lines: 29-43





