%% Create an instance of the class
clc;
L0              = 1;
Radius          = 0.01;
MassDensity     = 1062;
YoungModulus    = 6.66e5;
PoissonRatio    = 0.5;
MaterialDamping = 0.1;
NGaussPoints    = 20;
n               = 1;
StrainBasis     = @(s) [0;-1/L0;0;0;0;0];
Parameters      = [L0, Radius, MassDensity, YoungModulus, PoissonRatio, MaterialDamping]';
b1              = GVSBodyPlanar([Parameters; NGaussPoints]);
j1              = GVSJointPlanar([Parameters; NGaussPoints]);
%
%Build an equivalent body for comparison
b2              = PCCPlanar_Body_no_elongation(Parameters);
j2              = PCCPlanar_Joint_no_elongation(Parameters(1));
%Test configurations
q_test          = 1;
dq_test         = 1;
ddq_test        = -1;

%% Test all the methods
class_test;

%% Test dynamic methods
methods_test;


