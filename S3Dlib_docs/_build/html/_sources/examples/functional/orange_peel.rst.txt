.. _orange_peel:

.. role::  raw-html(raw)
    :format: html

*******************************************
Surface Texture
*******************************************

Surface textures may be generated using random color or geometry variations.
These effects can be accentuated using shading and highlighting.

Geometry
=========================

.. image:: images/orange_peel.png
   :class: sphx-glr-single-img

This first example combines several operations to control the final surface visualization. 

   1. randomly generated texture is applied.
   2. geometric functional operations are applied twice, in sequence.
   3. surface shading and highlighting are both applied.
   4. shading and highlighting are applied from two separate illumination directions.

A random generator was used to *slightly* change the coordinates of the unit radius
cylinder.  The twist function is similar to those in previous examples; however
for this case, the radial coordinate is used in the highlighted line instead of
unity.  Shading was then applied using the *direction = (1,1,1)*.  Because
the hilite direction is changed from default, the surface appears illuminated
from two light sources.  Texture enhances the effect.

For this example using a uniform surface color, the number of twists is
required to be even (2).  The reason being that shading is used to visualize the geometry.
In previous :ref:`twistribbon` examples, for an odd number of twists,
a mirrored cmap was needed so that the two ends at 0 
and :raw-html:`2&pi;` appear joined. Only when the twists are even are 
the surface normals continuous, which shading is based on. 


.. literalinclude:: source/ex_orange_peel.py
   :language: python
   :emphasize-lines: 20


Color Mapping
=========================

For this example, the random function, *randfunc*, is used to randomly select a color from the
*Oranges* Matplotlib colormap.

.. image:: images/orange_peel_cmap1.png
   :class: sphx-glr-single-img

The magnitude of the *sigma* parameter in the random function has no effect since the colormap
is normalized.  Each face color is randomly chosen and the resolution of the color fluctuations
are only dependent on the resolution set in the surface object.  The only change from
the previous script is higlighted.

.. literalinclude:: source/ex_orange_peel_cmap.py
   :language: python
   :lines: 22-28
   :emphasize-lines: 4


Combined Effects
=========================

In this example, the colormaping is based on the radial position of the randomized
surface.

.. image:: images/orange_peel_cmap2.png
   :class: sphx-glr-single-img

As in the first example, the initial surface is slightly changed by randomizing the
position of the radial coordinate.  This coordinate is then color mapped following
the randomization.  Note that the color mapping function now uses the radial 
coordinate instead of the random function, as was the case in the last example.

Since the colormapping uses the face coordinates which are computed from the three randomized
vertex coordinates, the effect is to smooth out the colorized texture.

.. literalinclude:: source/ex_orange_peel_cmap.py
   :language: python
   :lines: 30-37
   :emphasize-lines: 4,5






