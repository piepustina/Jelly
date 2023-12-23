PCC2D
----------

.. automodule:: classes

This class represents a planar slender soft body under the Piecewise Constant Curvature (PCC) hypothesis without elongation. The strain is modeled as

.. math:: 
     \boldsymbol{\xi}(q, s) = \boldsymbol{\Phi}(s)q + \boldsymbol{\xi}_{0} = \frac{1}{L_0} \left(\begin{array}{c} 1 \\ 0 \\ 0 \\ 0 \\ 0 \\ 1 \end{array}\right)q + \left(\begin{array}{c} 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 1  \end{array}\right),

where :math:`q \in \mathbb{R}` is the configuration variable and :math:`s \in [0, L_{0}]` the curvilinear abscissa, with :math:`L_{0}` the body rest length.

.. autoclass:: PCC2D
    :show-inheritance:
    :members: