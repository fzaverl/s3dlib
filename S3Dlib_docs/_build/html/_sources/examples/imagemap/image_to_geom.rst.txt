.. _image_to_geom:


**************************************
Image Geometric Mapping
**************************************

.. image:: images/image_to_geom.png
   :class: sphx-glr-single-img


The image pixel data is used to map color values to geometry.  For the left image, the value color component was
used (the default component).  For the right image, the hue color component was used by setting *cref* to 'h'.
Using the *hzero* parameter, the start value is blue (0.7) and increases in a negative
direction going from blue to cyan, green, red, etc.

.. image:: images/image_to_geom_sc.png
   :class: sphx-glr-single-img


.. literalinclude:: source/ex_image_to_geom.py
   :language: python

