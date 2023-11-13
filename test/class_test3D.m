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
n               = 3;
N_B             = 1;
Parameters      = [L0, BaseRadius, TipRadius, MassDensity, YoungModulus, PoissonRatio, MaterialDamping]';
b1 = PCC3D([Parameters; NGaussPoints]);
j1 = FixedJoint();
B1 = cell(N_B, 1);
J1 = cell(N_B, 1);
for i = 1:N_B
    B1{i} = PCC3D([Parameters; NGaussPoints]);
    J1{i} = FixedJoint();
end

%Build an equivalent body for comparison
Parameters      = [Parameters(1:2); Parameters(4:end)];
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


