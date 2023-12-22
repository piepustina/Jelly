VGVSBody
----------

.. automodule:: classes

This abstract class represents a soft body modeled under the Volumetric Geometric Variable Strain (VGVS) approach of :cite:p:`pustina2023univeral`.
The class implements all the methods required by a :class:`Body <Body>`. Concrete subclasses must define a strain basis :math:`\boldsymbol{\Phi}(\boldsymbol{q}, s) \in \mathbb{R}^{6 \times n}` such that 

.. math:: 
     \boldsymbol{\xi}(\boldsymbol{q}, s_{1}) = \boldsymbol{\Phi}(\boldsymbol{q}, s_{1})\boldsymbol{q} + \boldsymbol{\xi}_{0}

where :math:`\boldsymbol{\xi} \in \mathbb{R}^{6 \times 1}` is the strain, being :math:`\boldsymbol{\xi}_{0}` its stress-free value, :math:`\boldsymbol{q} \in \mathbb{R}^{n}` the vector of configuration variables and :math:`s_{1} \in [0, L_{0}]` the curvilinear abscissa, with :math:`L_{0}` the body rest length.

Concrete subclasses can also override the radius function to model the body of the radius as follows

.. math:: 
     R(\boldsymbol{q}, s_{1}, s_{2}) = R_{0} + \delta R(\boldsymbol{q}, s_{1}, s_{2}),

where :math:`R_{0}` is the stress-free radius and :math:`\delta R(\boldsymbol{q}, s_{1}, s_{2})` models the radius change.

See :class:`VPCC2D <VPCC2D>` for a possible implementation. 

.. autoclass:: VGVSBody
    :show-inheritance:
    :members:

.. toctree::
   
   VGVS/vpcc2d