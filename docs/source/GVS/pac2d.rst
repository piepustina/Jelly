PAC2D
---------------------------------------------------------

.. automodule:: classes

This class represents a slender planar soft body under the Piecewise Affine Curvature (PAC) hypothesis without elongation. The strain is modeled as

.. math:: 
     \boldsymbol{\xi}(\boldsymbol{q}, s) = \boldsymbol{\Phi}(s)\boldsymbol{q} + \boldsymbol{\xi}_{0} = \frac{1}{L_0} \left(\begin{array}{cc} 0 & 0 \\ -1 & -s \\ 0 & 0 \\ 0 & 0 \\ 0 & 0 \\ 0 & 0  \end{array}\right)\boldsymbol{q} + \left(\begin{array}{c} 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 1  \end{array}\right),

where :math:`\boldsymbol{q} \in \mathbb{R}^{2 \times 1}` is the vector of configuration variables and :math:`s \in [0, L_{0}]` the curvilinear abscissa, with :math:`L_{0}` the body rest length.

.. autoclass:: PAC2D
    :show-inheritance:
    :members: