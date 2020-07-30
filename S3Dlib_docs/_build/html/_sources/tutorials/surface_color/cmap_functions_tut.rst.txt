.. _cmap_functions_tut:

.. role::  raw-html(raw)
    :format: html


********************************
Color Mapping Using Functions
********************************


Surface color can be defined from function values that are mapped to colors using colormaps by calling the 
*surface* object method::

     surface.map_cmap_from_op( op, camp )


The *op* argument is a function with coordinate arguments and returning a single
value.  The  *cmap* argument is a color map or registered color map name.  Returned
function values at face coordinates are normalized and colors mapped from the colormap.

For this tutorial, the following two example functions are used to demonstrate this
method:

.. raw:: html

   f(x,y) = &half;sin( 6r )  <i>where</i> r = ( x<sup>2</sup> + y<sup>2</sup> ) <sup>&half;</sup>
   <br/>
   g(x,y) = sin( -xy )


First, consider the 2D visualization of the function f(x,y) using a simple Matplotlib script.

.. image:: images/col_func0.png
   :class: sphx-glr-single-img

The script, in the following, is organized in the three part structure to  
parallel the construction which will be further used for 3D visualizations.  Without a colorbar,
there is no indication as to the functional range for the applied domain.


.. literalinclude:: source/tut_col_func0.py
   :language: python


For a 3D visualization of the function, first consider applying a colormap based on one coordinate 
position of the surface, in this case the z-axis.  This is the 3D surface color mapping method
used in the numerous Matplotlib 
`3D plotting Gallery <https://matplotlib.org/3.1.1/gallery/index.html#d-plotting>`_
examples.

The following 3D plots of this function only differ by the projected views of the surface.  The default
Matplotlib view is used for the left plot.  The right plot used a view directly normal to
the xy plane (looking directly down from the top) .

For the left plot, the functional values are represented by the elevation of the surface
relative to the z-axis.  The colormap is actually redundant, it contributes very little
additional information. Surface shading would be sufficient for the 3D plot since the Z-axis
provides the functional values over the domain.

.. image:: images/col_func1.png
   :class: sphx-glr-single-img


When using the geometry to map colors, a mapping function can be defined which takes 3xN array
and returns a single N array value as::

    def func(xyz) :
        x,y,z = xyz
        # define v = f(x,y,z)
        return v

A lambda function may be convenient to return a single N value from a function.
For example::

    def otherfunc(xyz) :
         x,y,z = xyz
         # define X,Y,Z = f(x,y,z)
         return X,Y,Z

    func = lambda xyz : otherfunc(xyz)[2]

When the return value is simply based on the surface geometry coordinates, the lambda function
is simplified.  For example, when only the z-coordinate is used to map color, ie. z=xyz[2]::

     func = lambda xyz : xyz[2]

In the following script, the highlighted lines show the object method call used to
map the colors.  Instead of defining a separate function for the color mapping, lambda
functions were used as the mapping function argument.  The two different lambda functions
result in the identical values.  The first calls the geometry function and only returns
the z coordinate value.
The second lambda function just returns the z coordinate from the surface coordinates
which has already been mapped.  In both cases, the object color map is assigned 
during instantiation and then used in the functional mapping.


.. note::
   If the color map argument is not provided in the method call, the surface
   object color map is used.  If no color map is assigned during object instantiation,
   the Matplotlib default colormap of 'viridis' is used.


.. literalinclude:: source/tut_col_func1.py
   :language: python
   :emphasize-lines: 17,21,35,36


A core concept using S3Dlib is that surface geometry and surface color are controlled separately.  Whereas geometry
can be used to define color using a colormap, versatility of S3Dlib is using this separation.  This is
particularly useful as shown in the examples of :ref:`complex` or :ref:`sph_harm`.  In addition, surface
shape defined by parametric values were indicated by color in the  :ref:`anim_cat2heli` example.
In the examples section, surface shading is primarily only applied to enhance the visualization of the geometry,
whereas color is generally used for additional functional relationships.

.. image:: images/col_func2.png
   :class: sphx-glr-single-img

This concept is demonstrated above using both the f and g functions, but separately for
color and geometry.  In the two pair of highlighted lines in the script below, only the functional calls are reversed.
A moderated amount of shading was applied.  This was needed for the left plot to enhance the
geometry since only a slight variation in the g function occurs near the origin.  Further
annotation of these plots would be to include colorbar values.

.. literalinclude:: source/tut_col_func2.py
   :language: python
   :emphasize-lines: 21,22,25,26


