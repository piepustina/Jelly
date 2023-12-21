GVSBody
----------

.. automodule:: classes

This abstract class represents a slender soft body modeled under the Geometric Variable Strain (GVS) approach of :cite:p:`boyer2020dynamics`.
The class implements all the methods required by a :class:`Body <Body>`. Concrete subclasses must define a strain basis :math:`\boldsymbol{\Phi}(\boldsymbol{q}, s) \in \mathbb{R}^{6 \times n}` such that 

.. math:: 
     \boldsymbol{\xi}(\boldsymbol{q}, s) = \boldsymbol{\Phi}(\boldsymbol{q}, s)\boldsymbol{q} + \boldsymbol{\xi}_{0}

where :math:`\boldsymbol{\xi} \in \mathbb{R}^{6 \times 1}` is the strain, being :math:`\boldsymbol{\xi}_{0}` its stress-free value, :math:`\boldsymbol{q} \in \mathbb{R}^{n}` the vector of configuration variables and :math:`s \in [0, L_{0}]` the curvilinear abscissa, with :math:`L_{0}` the body rest length.

See :class:`PCC2D <PCC2D>` or :class:`PCC3D <PCC3D>` for how the strain basis can be implemented. 

.. autoclass:: GVSBody
    :show-inheritance:
    :members:

.. toctree::
   
   GVS/pcc2d
   GVS/pcc2dElongation
   GVS/pcc3d
   GVS/pac2d
   GVS/pac3d
   GVS/pgc3d