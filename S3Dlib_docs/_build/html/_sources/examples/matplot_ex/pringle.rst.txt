.. _pringle:


***************************
Color Mapping Normals
***************************

For this example, the surface geometry is similar to that used in the Matplotlib 
`Triangular 3D surfaces <https://matplotlib.org/3.1.1/gallery/mplot3d/trisurf3d.html#sphx-glr-gallery-mplot3d-trisurf3d-py>`_
example.  Note that the coordinates
are normalized, requiring scaling in the function definition.


.. image:: images/pringle.png
   :class: sphx-glr-single-img


.. literalinclude:: source/ex_pringle.py
   :language: python
 

To produce the similar colored view as the Matplolib plot, use shading instead of color mapping normals.
This is done by replacing the second section of script by::

    # 2. Setup and map surfaces .........................................
    
    saddle = s3d.PolarSurface(3, basetype='hex', linewidth=0.1)
    saddle.map_geom_from_op( pringle ).shade(0.3,direction=[-1,-1,1])

which produces the following plot.  In this case, a color is not defined in the constructor so the
blue default color is used for the surface color.

.. image:: images/pringle2.png
   :class: sphx-glr-single-img

