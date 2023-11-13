PCC2D
----------

.. automodule:: classes

Class representing a slender soft body modeled under the piecewise constant curvature hypothesis without elongation. The strain is modeled as

.. math:: 
     \xi(q, s) = \left(\begin{array}{c} 0 \\ \displaystyle -\frac{q}{L} \\ 0 \\ 0 \\ 0 \\ 1  \end{array}\right),

where :math:`L` is the rest length of the body.

.. autoclass:: PCC2D
    :show-inheritance:
    :members: