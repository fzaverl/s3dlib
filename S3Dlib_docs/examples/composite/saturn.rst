.. _saturn:


****************************
Alpha Channel Adjustments
****************************

.. image:: images/saturn.png
   :class: sphx-glr-single-img


This example demonstrates that the 'alpha' argument in the *set_surface_alpha* method
is a multiplier of the alpha color component and **not** a constant setting for the
alpha component.  It is set constant only if the 'constant' argument is set *True*
(if all pixels have the same alpha, the assignment of the constant is not needed).

The image, which sets the color of the rings, uses an alpha of 0 for pixels that
separate individual rings. Using *set_surface_alpha*, these pixels remain 0 when the
transparency for the fully opaque colors is reduced to 0.1.

.. literalinclude:: source/ex_saturn.py
   :language: python

