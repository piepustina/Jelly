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
N_B             = 2;
Parameters      = [L0, BaseRadius, TipRadius, MassDensity, YoungModulus, PoissonRatio, MaterialDamping]';
b1 = PCC3D([Parameters; NGaussPoints]);
j1 = FixedJoint();
B1 = cell(N_B, 1);
J1 = cell(N_B, 1);
for i = 1:N_B
    B1{i} = PCC3D([Parameters; NGaussPoints]);
    %B1{i} = PCC2D([Parameters; NGaussPoints]);
    J1{i} = FixedJoint();
end

%% Test the direct and differential kinematics with respect to the bodytree class
r1 = SoftRobot(J1, B1, {});
r2 = BodyTree(J1, B1);

%%
q_test  = ones(r1.n, 1);
dq_test = ones(r1.n, 1);

% Transformation matries
T1 = r1.DirectKinematics(q_test, L0*(1:N_B))

%T2 = r2.DirectKinematics(q_test)

% Body Jacobian
J1 = r1.BodyJacobian(q_test, L0*(1:N_B))

%J2 = r2.BodyJacobian(q_test)

%% Test the inverse kinematics with a random configuration
q_min  = -2;
q_max  = 2;
disp("The configuration is ")
q_test = (q_max-q_min).*rand(r1.n,1) + q_min
% Direct kinematics
T = r1.DirectKinematics(q_test, L0*(1:N_B));

% Inverse kinematics
disp("The inverse kinematics result is ")
tic
[q_ik, converged] = r1.InverseKinematics(T, L0*(1:N_B), zeros(r1.n, 1), 5)
toc

disp("The error norm is " + norm(q_test - q_ik))

%% Test also the simulink model
open("simulink_soft_ik.slx")
