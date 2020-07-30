.. _cube:


***********************
Base Class Surface
***********************

.. image:: images/cube.png
   :class: sphx-glr-single-img


Most examples use object classes derived from the Surface3DCollection base class.
This example demonstrates using 'raw' vertex and face arrays to instantiate a base class object,
which as a result, uses xyz as native coordinates.  When using this base class, all faces must
have the same number of vertices (3,4 or 5).  Faces may be further subdivided into triangles using the
*triangulate()* method, as shown in :ref:`intermediate` example.



.. literalinclude:: source/ex_cube.py
   :language: python

