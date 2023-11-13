GVSBody
----------

.. automodule:: classes

This abstract class represents a slender soft body modeled under the Geometric Variable Strain (GVS) hypothesis. 
The class implements all the methods required by a :class:`Body <Body>`, but requires the definition of a strain function :math:`\xi(q, s)`, where :math:`q` is the vector of configuration variables
and :math:`s \in [0, L]` the arc length, with :math:`L` the rest length of the body. 

See :class:`PCC2D <PCC2D>` or :class:`PCC3D <PCC3D>` for possible implementations. 

.. autoclass:: GVSBody
    :show-inheritance:
    :members:

.. toctree::

   PCC/pcc3d
   PCC/pcc2d