%Inizialize the Soft Robotics toolbox
initToolbox;

%%
%Compute all the required terms of each body
%PCC_planar_slender_model_no_elongation;
%PCC_planar;

%%

n_B = 2;
L_0_s = 0.3*ones(n_B, 1);
L_0_s = ones(n_B, 1);

R     = 0.01*ones(n_B, 1);
rho_s   = 1062*ones(n_B, 1);
E     = 6.66e5*ones(n_B, 1);
Poi   = 0.2*ones(n_B, 1);
Eta   = 1e4*ones(n_B, 1);

%Define the joints of the robot
Joints = cell(n_B, 1);
Bodies = cell(n_B, 1);

for i = 1:n_B
    Parameters = [L_0_s(i), R(i), rho_s(i), E(i), Poi(i), Eta(i)];
    Joints{i} = PCCPlanar_Joint(Parameters(1));
    Bodies{i} = PCCPlanar_Body(Parameters);
end

%Construct a BodyTree of the robot
bt_2D = BodyTree(Joints, Bodies, {});

%Store the total number of DOFs
n = bt_2D.n;

%Generate a random reference in between "a" and "b"
rng(0,'twister');
a = -1;
b = 1;
ref = (b-a).*rand(n, 1) + a;

%% Try to compute the forward and inverse dynamics

bt_2D.ForwardDynamics(zeros(n, 1), ones(n, 1), 0, 'double')

%%

bt_2D.InverseDynamics(ones(n, 1), ones(n, 1), zeros(n, 1), 'double')



