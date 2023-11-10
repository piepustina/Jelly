%% Create an instance of the class
clear; clc;
L0              = 0.5;
Radius          = 0.01;
MassDensity     = 1062;
YoungModulus    = 6.66e5;
PoissonRatio    = 0.5;
MaterialDamping = 0.1;
NGaussPoints    = 10;
n               = 3;
N_B             = 2;
StrainBasis     = @(s) [   0,   -1/L0, 0;
                        1/L0,       0, 0;
                           0,       0, 0;
                           0,       0, 0;
                           0,       0, 0;
                           0,       0, 1/L0];
Parameters      = [L0, Radius, MassDensity, YoungModulus, PoissonRatio, MaterialDamping]';
b1 = GVSBody3D([Parameters; NGaussPoints]);
j1 = FixedJoint();
B1 = cell(N_B, 1);
J1 = cell(N_B, 1);
for i = 1:N_B
    B1{i} = GVSBody3D([Parameters; NGaussPoints]);
    J1{i} = FixedJoint();
end

%Build an equivalent body for comparison
b2 = PCCCylindricalDelta_Body(Parameters);
j2 = PCCCylindricalDelta_Joint(Parameters(1));
B2 = cell(N_B, 1);
J2 = cell(N_B, 1);
for i = 1:N_B
    B2{i} = PCCCylindricalDelta_Body(Parameters);
    J2{i} = PCCCylindricalDelta_Joint(Parameters(1));
end
%Test configurations
q_test      = 1*[1; -1; 0.2];
dq_test     = 0*[1; 0.3; 1];
ddq_test    = 0*[-1; 1; 0.5];

%% Test all the methods
class_test;

%% Test dynamic methods
methods_test;


