PGC3D
----------

.. automodule:: classes

This class represents a slender soft body under the Piecewise Gaussian Curvature (PGC) hypothesis with elongation. The strain is modeled as

.. math:: 
     \boldsymbol{\xi}(\boldsymbol{q}, s) = \boldsymbol{\Phi}(s)\boldsymbol{q} + \boldsymbol{\xi}_{0} = \frac{1}{L_0} \left(\begin{array}{ccccc} 0 & 0 & -1 & -e^{-(s_{i}-L_{0_{i}}/2)^2} & 0 \\
        1 & e^{-(s_{i}-L_{0_{i}}/2)^2} & 0 & 0 & 0\\
        0 & 0 & 0 & 0 & 1  \end{array}\right)\boldsymbol{q} + \left(\begin{array}{c} 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 1  \end{array}\right),

where :math:`\boldsymbol{q} \in \mathbb{R}^{5 \times 1}` is the vector of configuration variables and :math:`s \in [0, L_{0}]` the curvilinear abscissa, with :math:`L_{0}` the body rest length.

The class is used in Simulation 2 of :cite:p:`pustina2023univeral`. 

.. autoclass:: PGC3D
    :show-inheritance:
    :members: