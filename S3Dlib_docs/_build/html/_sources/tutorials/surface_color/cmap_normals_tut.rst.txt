.. _cmap_normals_tut:

****************************************
Color Mapped Surface Normals
****************************************


Surface color can be defined based on the surface normals that are mapped to colors by
calling the *surface* object method::

   surface.map_cmap_from_normals( cmap, direction )

The *cmap* argument is a color map or registered color map name. The *direction* argument
is a 3D array in x,y,z coordinates, which is the direction of incident light.
If no arguments are provided, as shown in the following plot, the default values
for cmap and direction are *viridis*, the Matplot default, and [1,0,1], respectively.

.. image:: images/col_norm0.png
   :class: sphx-glr-single-img



The only difference in the :ref:`hello-1` script from the script given in the *Hello World* tutorial 
is that normal color mapping is used instead of shading, as shown in the highlighted 
line.



.. literalinclude:: source/tut_col_norm0.py
   :language: python
   :emphasize-lines: 17


By assigning a *cmap* value, various coloring effects can be created.  For example,
using::

   surface.map_cmap_from_normals( 'magma' )

The following plot was created.

.. image:: images/col_norm1.png
   :class: sphx-glr-single-img


More details are provided in the :ref:`color_mapping_normals` guide.
 
Shading is the only color mapping operation that is additive.  Shading can be applied
to any surface.

.. image:: images/col_norm2.png
   :class: sphx-glr-single-img

.. literalinclude:: source/tut_col_norm1.py
   :language: python
   :emphasize-lines: 16,23

