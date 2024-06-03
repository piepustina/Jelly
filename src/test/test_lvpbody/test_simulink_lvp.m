%% Create an instance of the class
clear; clc;

N_B = 4;

% Load a test mesh
%femodel = femodel(Geometry="./test/meshes/Cylinder_coarse.stl");
%femodel = femodel(Geometry="./test/meshes/Diamond.stl");
%femodel = femodel(Geometry="./test/meshes/Diamond_high_res.stl");
femodel = femodel(Geometry="./test/meshes/Diamond_low_res.stl");



model = generateMesh(femodel, "Hmax", 50, "Hmin", 50);
%model = generateMesh(femodel, "Hmax", 10, "Hmin", 9);
%model = generateMesh(femodel);

Nodes    = model.Geometry.Mesh.Nodes./1000;% Scale from [mm] to [m]. The mesh should be already with [m] units!
Elements = model.Geometry.Mesh.Elements;

% Parameters of the body
L0              = max(Nodes(3, :));
nGauss          = 2;
MassDensity     = 960;
YoungModulus    = 5e5;
PoissonRatio    = 0.4;
DampingFactor   = 0.05;

%% Build the bodies
B1 = cell(N_B, 1);
J1 = cell(N_B, 1);
for i = 1:N_B
    Primitives = {PCStretchCompressionPrimitive(L0); ...
                  PCSTBPrimitive(L0)};
    B1{i} = LVPBody(Nodes, Elements, Primitives, nGauss, MassDensity, YoungModulus, PoissonRatio, DampingFactor);
    J1{i} = FixedJoint();
end

B1{end+1} = RigidBody([10; zeros(9, 1)]);
J1{end+1} = FixedJoint();


%% Build the robot
%r1 = BodyTree(J1, B1);
r1 = LVPBodyTree(J1, B1);
%r1 = LVPBodyTreeJAct(J1, B1);
%r1.g = [0; -9.81; 0];
%r1.g = [-9.81; 0; 0];
r1.g = [0; 0; 9.81];

%% Test the parallelization
a = -10;
b =  10;
q = (b-a).*rand(r1.n, 1) + a;

tic,
[M1, t1 ] = r1.UnifiedInverseDynamics(q, ones(r1.n, 1), ones(r1.n, 1));
toc
tic
[M2, t2 ] = r1.ParUnifiedInverseDynamics(q, ones(r1.n, 1), ones(r1.n, 1));
toc

%% Solve for the equilibrium
close all;
figure; hold on; grid on; view(3)
light("Position", [-0.1, -0.1, 0.1])
lighting gouraud

[q_ss, ~] = r1.EquilibriumConfiguration(zeros(r1.n, 1), [200; 0; 0; 0])
r1.plot(q_ss, "LineStyle", "-", "FaceAlpha", 1);
%r1.plot([0; 0; 0; 0; pi/8; 0; 0; 0; 0; 0; 0; 0], "LineStyle", "-", "FaceAlpha", 1);

xlabel("$x [m]$", "Interpreter", "latex", "FontSize", 14)
ylabel("$y [m]$", "Interpreter", "latex", "FontSize", 14)
zlabel("$z [m]$", "Interpreter", "latex", "FontSize", 14)

axis equal

%% Open the simulink system
open("simulink_test_bt_lvp.slx")

%% Plot the final configuration
close all;
figure; hold on; grid on; view(3)
light("Position", [-0.1, -0.1, 0.1])
lighting gouraud

r1.plot(squeeze(out.q.Data(:, :, end)), "LineStyle", "-", "FaceAlpha", 1);

xlabel("$x [m]$", "Interpreter", "latex", "FontSize", 14)
ylabel("$y [m]$", "Interpreter", "latex", "FontSize", 14)
zlabel("$z [m]$", "Interpreter", "latex", "FontSize", 14)

axis equal