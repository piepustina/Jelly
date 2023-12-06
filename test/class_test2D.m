%% Create an instance of the class
clear; clc;
L0              = 1;
BaseRadius      = 0.01;
TipRadius       = 0.01;
MassDensity     = 1062;
YoungModulus    = 6.66e5;
PoissonRatio    = 0.5;
MaterialDamping = 0.1;
NGaussPoints    = 10;
n               = 1;
N_B             = 1;
Parameters      = [L0, BaseRadius, TipRadius, MassDensity, YoungModulus, PoissonRatio, MaterialDamping]';
b1              = PCC2D([Parameters; NGaussPoints]);
j1              = FixedJoint();
B1 = cell(N_B, 1);
J1 = cell(N_B, 1);
for i = 1:N_B
    B1{i} = PCC2D([Parameters; NGaussPoints]);
    J1{i} = FixedJoint();
    %J1{i} = RotationalJoint([0; 0; 0; 0]);
end
%
%Build an equivalent body for comparison
Parameters      = [Parameters(1:2); Parameters(4:end)];
b2              = PCCPlanar_Body_no_elongation(Parameters);
j2              = PCCPlanar_Joint_no_elongation(Parameters(1));
B2 = cell(N_B, 1);
J2 = cell(N_B, 1);
for i = 1:N_B
    B2{i} = PCCPlanar_Body_no_elongation(Parameters);
    J2{i} = PCCPlanar_Joint_no_elongation(Parameters(1));
end

%Test configurations
q_test          = 1;
dq_test         = 1;
ddq_test        = -1;

%% Test all the methods
class_test;

%% Test dynamic methods
methods_test;


