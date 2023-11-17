%% Create an instance of the class
clear; clc;
L0              = 0.5;
BaseRadius      = 0.01;
TipRadius       = 0.01;
MassDensity     = 1062;
YoungModulus    = 6.66e5;
PoissonRatio    = 0.5;
MaterialDamping = 0.1;
NGaussPoints    = 10;
Parameters      = [L0, BaseRadius, TipRadius, MassDensity, YoungModulus, PoissonRatio, MaterialDamping]';
B = cell(3, 1);
J = cell(3, 1);

B{1} = PCC3D([Parameters; NGaussPoints]);
J{1} = FixedJoint();
B{3} = PCC3D([Parameters; NGaussPoints]);
J{3} = FixedJoint();

%Rigid Body parameters
m               = 1;
I               = eye(3);
pcom            = [1; 0; -1];
Parameters  = [m; pcom; I(1, 1); I(2, 2); I(3, 3); I(1, 2); I(1, 3); I(2, 3)];
B{2} = RigidBody(Parameters);
J{2} = FixedJoint();


%% Create the tree
r1           = BodyTree(J, B);

q_test      = [1; 1; 0.1; -1; 1; 0];
dq_test     = zeros(6, 1);



