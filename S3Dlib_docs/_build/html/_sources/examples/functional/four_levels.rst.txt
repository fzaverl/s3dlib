.. _four_levels:


****************************
Segmented Cmap Operation
****************************

.. image:: images/four_levels.png
   :class: sphx-glr-single-img

This example demonstrates two methods of appling identical surface coloring.  The only difference is the
resulting colorbar reference.

Colors were selected from a Matplotlib sequential colormap, 'plasma'.  This provided selecting colors which
are visualy different in lightness, L*, which are perceived differently using gray scale printing.

.. literalinclude:: source/ex_four_levels.py
   :language: python
   :emphasize-lines: 32

