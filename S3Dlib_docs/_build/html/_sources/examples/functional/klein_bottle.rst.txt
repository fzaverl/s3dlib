.. _klein_bottle:


*******************************************
Klein Bottle, Spherical to XYZ
*******************************************

.. image:: images/klein_bottle.png
   :class: sphx-glr-single-img

The main point of this example is that

*The development of the functional definition of a surface is the hard part, visualizing the
three dimensional surface is fairly easy using S3Dlib with Matplotlib.*

As seen in the below script, this surface was constructed using a SphericalSurface object.
Alternatively, a PlanarSurface object could be used as demonstrated in
the :ref:`planarKlein` guide plot.

A detailed description of a Klein Bottle is found in
`Wikipedia <https://en.wikipedia.org/wiki/Klein_bottle#Bottle_shape>`_
where the functional definition is located.

.. literalinclude:: source/ex_klein_bottle.py
   :language: python

