Basic usage
===========

.. _installation:

Installation
------------

Download the source code from its GitHub page

.. code-block:: console

   (.venv) $ pip install lumache

and place it in your working directory. 

.. _initialization:

Initialization
--------------

Add the library to the MATLAB path by running the ``initToolbox.m`` script. 

.. code-block:: matlab
   :linenos:

   initToolbox.m

.. note::
   The script has to be executed every time a new instance of MATLAB is started. 

Introduction to Jelly
---------------------
Suppose we want to model a 3D continuum soft robot of two bodies modeled under the Piecewise Constant Curvature (PCC) hypothesis.
Each body has rest length :math:`L_{0} = 0.3 [\mathrm{m}]`, base and tip radius :math:`R_{\mathrm{base}} = 0.01 [\mathrm{m}]` and :math:`R_{\mathrm{tip}} = 0.01 [\mathrm{m}]`, respectively, 
and mass density :math:`\rho = 1062 [\mathrm{kg/m^{3}}]`. 
Furthermore, :math:`E = 6.66 [\mathrm{MPa}]` is the Young modulus, :math:`\nu = 0.5` the Poisson ratio and :math:`\eta = 0.1 [\mathrm{s}]` the material damping coefficient.
Finally, use :math:`N_{\mathrm{Gauss}} = 10` Gaussian for the numerical computation of the kinematics and the integrals.

To create such robot two main ingredients are needed:

#. Two joints that connect the robot bodies;
#. Two 3D PCC bodies.

Create two fixed joints to connect the first body to the robot base and second body to the first.

.. code-block:: matlab
   :linenos:

   %Allocate the joints and store them in a cell array
   Joints = {FixedJoint(); FixedJoint()};

Now, create two 3D PCC bodies by using the :class:`PCC3D <PCC3D>` class.

.. code-block:: matlab
   :linenos:

   %Rest length
   L              = 0.3;
   %Base radius
   R_base         = 0.01;
   %Tip radius
   R_tip          = 0.01;
   %Mass density
   rho            = 1062;
   %Young modulus
   E              = 6.66*1e6;
   %Poisson ratio
   nu             = 0.5;
   %Material damping
   eta            = 0.1;
   %Number of Gaussian points
   N_Gauss        = 10;
   %Create a vector of parameters
   BodyParameters = [L; R_base; R_tip; rho; E; nu; eta; N_Gauss];
   %Allocate the bodies and store them in a cell array
   Bodies         = {PCC3D(BodyParameters); PCC3D(BodyParameters)};

In **Jelly**, a robot is modeled by the :class:`BodyTree <BodyTree>` class, which is a serial assembly of :math:`N` :class:`joints <Joint>` and :class:`bodies <Body>`. 
To build the robot create a :class:`BodyTree <BodyTree>` instance.

.. code-block:: matlab
   :linenos:

   Robot = BodyTree(Joints, Bodies);

Now, use the ``Robot`` object to compute some dynamic terms.  

.. code-block:: matlab
   :linenos:

   %Configuration variables and time derivatives
   q   = [0; -pi; 0; pi/3; pi/4; -0.01];
   dq  = zeros(Robot.n, 1);
   ddq = zeros(Robot.n, 1);
   %Generalized actuation force
   tau = ones(Robot.n, 1);

   %Evaluate the inverse dynamics
   Robot.InverseDynamics(q, dq, ddq)

   %Evaluate the forward dynamics
   Robot.ForwardDynamics(q, dq, tau)

   %Evaluate the mass matrix
   Robot.MassMatrix(q)

The ``Robot`` can also be used in Simulink for real-time control or simulate the dynamics. 
Create a blank Simulink model and open it. In the library browser, 