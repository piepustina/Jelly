BodyTree
----------

.. automodule:: classes

This class represents a serial modular mechanical system :cite:p:`pustina2023univeral`. 
A :class:`BodyTree <BodyTree>` consists of :math:`N` :class:`joints <Joint>` and :class:`bodies <Body>` serially interconneted. 

The dynamic terms are computed through the Generalized Inverse Dynamics (GID) and Actuation Inverse Dynamics (AID) algorithms of :cite:p:`pustina2023univeral`.

.. autoclass:: BodyTree
    :show-inheritance:
    :members: