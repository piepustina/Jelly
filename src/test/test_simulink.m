%% Create an instance of the class
clear; clc;
L0              = 1;
BaseRadius      = 0.01;
TipRadius       = 0.01;
MassDensity     = 1062;
YoungModulus    = 6.66e5;
PoissonRatio    = 0.5;
MaterialDamping = 0.1;
NGaussPoints    = 5;
n               = 1;
N_B             = 3;
Parameters      = [L0, BaseRadius, TipRadius, MassDensity, YoungModulus, PoissonRatio, MaterialDamping]';
b1              = PCC2D([Parameters; NGaussPoints]);
j1              = FixedJoint();
B1 = cell(N_B, 1);
J1 = cell(N_B, 1);
for i = 1:N_B
    %B1{i} = PCC2D([Parameters; NGaussPoints]);
    B1{i} = PCC2DElongation([Parameters; NGaussPoints]);
    J1{i} = FixedJoint();
end

%Create the actuators
A2 = cell(N_B, 1);
for i = 1:N_B
   Parameters  = [(i-1)*L0; i*L0; NGaussPoints; BaseRadius/2; 0; 0];
   A2{i}       = ConstantDistanceActuator(Parameters);
end

A1 = cell(2, 1);
for i = 1:2
   if ~mod(i, 2)
       signA = -1;
   else
       signA = 1;
   end
   Parameters  = [0; N_B*L0; NGaussPoints; signA*BaseRadius/2; 0; 0];
   A1{i}       = ConstantDistanceActuator(Parameters);
end


%% Build the robot
%r1 = BodyTree(J1, B1);

r1 = SoftRobot(J1, B1, A1);

r2 = SoftRobot(J1, B1, A2);

%% Open the simulink system
open("simulink_test.slx")

