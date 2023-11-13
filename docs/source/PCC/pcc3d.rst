PCC3D
----------

.. automodule:: classes

Class representing a slender soft body modeled under the piecewise constant curvature hypothesis with elongation. 
The strain is modeled as

.. math:: 
     \xi(q, s) = \left(\begin{array}{c} \displaystyle \frac{q_{2}}{L} \\ \displaystyle -\frac{q_{1}}{L} \\ 0 \\ 0 \\ 0 \\ \displaystyle \frac{q_{3} + L}{L}  \end{array}\right),

where :math:`L` is the rest length of the body.

.. autoclass:: PCC3D
    :show-inheritance:
    :members: