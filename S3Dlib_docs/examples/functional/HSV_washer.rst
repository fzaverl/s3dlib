.. _hsv_washer:


**********************************************
Functional HSV Color Mapping
**********************************************

.. image:: images/hsv_washer.png
   :class: sphx-glr-single-img


This example is based on the Matplotlib function used in the
`cylindrical volumetric plot <https://matplotlib.org/gallery/mplot3d/voxels_torus.html#sphx-glr-gallery-mplot3d-voxels-torus-py>`_
example.


.. literalinclude:: source/ex_hsv_washer.py
   :language: python

A cylindrical surface is first mapped to a flat ring using the 'Ring()' function.  The figure at the below
left is a result of a direct call to this function, resulting in an intermediate surface.
Then, this surface only calls the 'warp()' function for geometric mapping.
This multiple mapping technique was used in the  :ref:`knot` example. 

The color is simply a result of using the rtz coordinates in hsv space.  The
variation of saturation in the radial direction is more apparent with a thicker surface as shown in the figure at
the lower right.  HSV space is more clearly illustrated in the :ref:`Lab_space` example and 
the :ref:`anim_hsv_cylinder` animation.

.. image:: images/hsv_washer2.png
   :class: sphx-glr-single-img


