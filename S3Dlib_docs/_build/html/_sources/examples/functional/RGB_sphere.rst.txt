.. _RGB_sphere:


**********************************************
Functional RGB Color Mapping
**********************************************

.. image:: images/RGB_sphere.png
   :class: sphx-glr-single-img


This example is based on the Matplotlib function used in the
`RGB volumetric plot <https://matplotlib.org/gallery/mplot3d/voxels_rgb.html#sphx-glr-gallery-mplot3d-voxels-rgb-py>`_
example.  RGB space is more clearly illustrated in the :ref:`Lab_space` example and the :ref:`anim_rgb_cube` animation.


.. literalinclude:: source/ex_RGB_sphere.py
   :language: python

.. note::
   Once transformed from rtp coordinates, the color is simply a result of using the XYZ coordinates in RGB space. 
   Since the sphere is initially expanded, the RGB colors are clipped to the range [0,1].


