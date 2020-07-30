.. _mp_wire:


**************************************
Datagrid Wireframe Plot
**************************************

Surface geometry is that taken from the matplotlib 
`3D wireframe plot <https://matplotlib.org/3.1.1/gallery/mplot3d/wire3d.html#sphx-glr-gallery-mplot3d-wire3d-py>`_
example.


.. image:: images/mp_wire.png
   :class: sphx-glr-single-img

When using a datagrid directly with a PlanarSurface object, the data is scaled along the
z-axis range [0,1].   Additional scaling is demonstrated in the :ref:`wireframe` example.



.. literalinclude:: source/ex_mp_wire.py
   :language: python
   :emphasize-lines: 29

