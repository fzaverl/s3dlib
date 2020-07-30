.. _triple:


**********************************************
Cmapped Normals, Shading and Highlighting
**********************************************

.. image:: images/triple.png
   :class: sphx-glr-single-img


This example is based on the Matplotlib function used in the
`Hillshading <https://matplotlib.org/3.1.1/gallery/specialty_plots/advanced_hillshading.html>`_
example.

.. note::
   In this example, the 3D geometry is used to represent the functional relationship
   and the **color is only used to visualize the geometry** .  This is in comparison
   to the *Matplotlib Hillshading* example where color is used for the functional
   value, hence the need of a colorbar.



.. literalinclude:: source/ex_triple.py
   :language: python

In a similar manner as the Hillshading example, two different functions are
used to apply color and geometry. 

.. image:: images/triple2.png
   :class: sphx-glr-single-img

The *wavefunc* is still used to apply color, but now the *surfripple* function is
used for the geometry.

.. literalinclude:: source/ex_triple2.py
   :language: python
   :lines: 7-29

