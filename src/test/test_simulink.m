%% Create an instance of the class
clear; clc;
L0              = 0.2;
BaseRadius      = 0.01;
TipRadius       = 0.01;
MassDensity     = 1062;
YoungModulus    = 6.66e5;
PoissonRatio    = 0.5;
MaterialDamping = 0.1;
NGaussPoints    = 10;
n               = 1;
N_B             = 3;
Parameters      = [L0, BaseRadius, TipRadius, MassDensity, YoungModulus, PoissonRatio, MaterialDamping]';
b1              = PCC2D([Parameters; NGaussPoints]);
j1              = FixedJoint();
B1 = cell(N_B, 1);
J1 = cell(N_B, 1);
for i = 1:N_B
    B1{i} = PCC2D([Parameters; NGaussPoints]);
    J1{i} = FixedJoint();
end


%% Build the robot
%r1 = BodyTree(J1, B1);

r1 = SoftRobot(J1, B1);

%% Open the simulink system
open("simulink_test.slx")

