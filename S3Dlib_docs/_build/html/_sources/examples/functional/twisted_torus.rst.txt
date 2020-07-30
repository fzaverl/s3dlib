.. _order_operation:

.. role::  raw-html(raw)
    :format: html

******************************************************
Order of Operation
******************************************************

.. image:: images/twisted_torus.png
   :class: sphx-glr-single-img

A basetype of 'squ_s' was used for this case 
because the twisted surface is rejoined at 
:raw-html:`&theta;` equal to 0 and :raw-html:`2&pi;`.
A custom color map, using the *cmap_utilities*, was used to emphasize the effect.


The color map was first applied to the cylinder.
Then the twisting geometric mapping was made.  The color map was applied first because the geometry
has the same 'shape' both before and after the twist operation is applied.  

.. literalinclude:: source/ex_twisted_torus.py
   :language: python
   :emphasize-lines: 32,33

If the highlighted lines in the above code are reversed, the visualization of
4 twists is lost, as seen below.

.. image:: images/twisted_torus2.png
   :class: sphx-glr-single-img


Notice the scale on the above two colorbars are not the same.  The top plot
colorbar is normalized from the top and bottom of the original cylindrical
surface prior to geometric mapping.  The bottom plot reflects the upper and
lower z-coordinate boundaries of the torus.  ( note the ratio = .45 in the
code which is reflected in this lower colorbar scale )