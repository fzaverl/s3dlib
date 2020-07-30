.. _hirez:


**************************************
High Resolution Example
**************************************

.. image:: images/hirez.png
   :class: sphx-glr-single-img

Mapping function definition was taken from the
`Mayavi <https://docs.enthought.com/mayavi/mayavi/mlab.html#id5>`_
demo example.  Notice in the script, the function only transforms the radial coordinate
since the native coordinates for a SphericalSurface object are already spherical coordinates.
Since the Mayavi example is symetric about the y-axis, a rotational transform was
applied to create a similar visual orientation.  Also, Numpy arrays are used in this functional
definition.

This S3Dlib example demonstrates the limitations using the predefined surface grids.  Here, the visualization
anomalies are along surface faces near the axis of symetry.  These anomalies could be eliminated with
an optimized grid specifically designed for this function.  The tradeoff is with development time which
can be reduced using S3Dlib surface objects.
A very minor effort is needed to create a custom colormap, apply the functional geometry, orient
the surface for viewing and adjust the shading/highlighting.  Then let Matplotlib do the heavy lifting of
putting the object on the screen.


.. note::

   This example image uses a high rez of 8.  Producing about 4 million colored faces takes a bit of execution
   time.  For the development of this script, a rez of 3 was used to determine the desired orientation and
   coloring.  Once set, the rez was set to 8 and the image computed.


.. literalinclude:: source/ex_hirez.py
   :language: python

