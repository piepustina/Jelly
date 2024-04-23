%% Load a test mesh
clc; clear; close all;
femodel = femodel(Geometry="./test/meshes/Cylinder_coarse.stl");

%model = generateMesh(femodel, "Hmax", 60, "Hmin", 60, "Hgrad", 2);
%model = generateMesh(femodel, "Hmax", 10, "Hmin", 9);
model = generateMesh(femodel);

% Get the nodes and elements
Nodes    = model.Geometry.Mesh.Nodes./1000;% Scale from [mm] to [m]. The mesh should be already with [m] units!
Elements = model.Geometry.Mesh.Elements;

%% Build the LVPBody
Primitives = {PCStretchCompressionPrimitive(0.3), ...
              PCTwistShearPrimitive(0.3), ...
              PCBendingPrimitive(0.3)};

B = LVPBody(Nodes, Elements, Primitives, 10, 1062, 5e4, 0.5, 0.05);

%% Compute the kinematics
% Configuration
q_min   = -2;
q_max   = 2;
q_test  = (q_max-q_min).*rand(B.n,1) + q_min;

q_test = [-1; 10; -0.1; 0.1; 2; 1];


% First order time derivative of the configuration
dq_min  = -100;
dq_max  =  100;
dq_test = (dq_max-dq_min).*rand(B.n,1) + dq_min;
% Second order time derivative of the configuration
ddq_min  = -5;
ddq_max  =  5;
ddq_test = (ddq_max-dq_min).*rand(B.n,1) + ddq_min;

tic
[xq, dxq, ddxq, Jq, Jx_ref, JJq] = B.UpdateKinematics(q_test, dq_test, ddq_test); 
toc;

% Check the error between the Jacobian and the time derivative of x
dxJ = squeeze(pagemtimes(Jq, dq_test));
rmse(dxq, dxJ, 2)

%%
close all
figure; hold on; grid on; view(3)
light("Position", [-0.1, -0.1, 0.1])
lighting gouraud
B.plot(q_test, "LineStyle", "-", "FaceAlpha", 1);

axis equal

%% Create the robot
Rb      = BodyTree({FixedJoint(); FixedJoint()}, {B; RigidBody([1; zeros(9, 1)])});
%Rb.g    = [0, -9.81, 0]';
Rb.g    = [0, 0, 9.81]';

q_eq    = Rb.EquilibriumConfiguration(zeros(6, 1) , zeros(6, 1));

%%
close all
figure; hold on; grid on; view(3)
light("Position", [-0.1, -0.1, 0.1])    
lighting gouraud
B.plot(q_eq, "LineStyle", "-", "FaceAlpha", 1);

axis equal


