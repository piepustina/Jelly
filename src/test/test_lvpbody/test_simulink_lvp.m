%% Create an instance of the class
clear; clc;

N_B = 1;

% Load a test mesh
%femodel = femodel(Geometry="./test/meshes/Cylinder_coarse.stl");
femodel = femodel(Geometry="./test/meshes/Diamond.stl");

%model = generateMesh(femodel, "Hmax", 60, "Hmin", 60, "Hgrad", 2);
model = generateMesh(femodel, "Hmax", 10, "Hmin", 9);
%model = generateMesh(femodel);

Nodes    = model.Geometry.Mesh.Nodes./1000;% Scale from [mm] to [m]. The mesh should be already with [m] units!
Elements = model.Geometry.Mesh.Elements;

% Parameters of the body
L0              = max(Nodes(3, :));
nGauss          = 2;
MassDensity     = 960;
YoungModulus    = 1e6;
PoissonRatio    = 0.5;
DampingFactor   = 0.05;

%% Build the bodies
B1 = cell(N_B, 1);
J1 = cell(N_B, 1);
for i = 1:N_B
    Primitives = {PCStretchCompressionPrimitive(L0), ...
                  PCTwistShearPrimitive(L0), ...
                  PCBendingPrimitive(L0)};
    B1{i} = LVPBody(Nodes, Elements, Primitives, nGauss, MassDensity, YoungModulus, PoissonRatio, DampingFactor);
    J1{i} = FixedJoint();
end

B1{end+1} = RigidBody([1; zeros(9, 1)]);
J1{end+1} = FixedJoint();


%% Build the robot
%r1 = BodyTree(J1, B1);
r1 = LVPBodyTree(J1, B1);
r1.g = [0; -9.81; 0];

%% Open the simulink system
open("simulink_test_bt_lvp.slx")

%% Plot the final configuration
close all;
figure; hold on; grid on; view(3)
light("Position", [-0.1, -0.1, 0.1])
lighting gouraud

B1{1}.plot(squeeze(out.q.Data(:, :, end)), "LineStyle", "-", "FaceAlpha", 1);

xlabel("$x [m]$", "Interpreter", "latex", "FontSize", 14)
ylabel("$y [m]$", "Interpreter", "latex", "FontSize", 14)
zlabel("$z [m]$", "Interpreter", "latex", "FontSize", 14)

axis equal