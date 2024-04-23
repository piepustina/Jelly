% Load a test mesh
clc; clear; close all;
femodel = femodel(Geometry="./test/meshes/Cylinder_coarse.stl");

model = generateMesh(femodel, "Hmax", 60, "Hmin", 60, "Hgrad", 2);
%model = generateMesh(femodel, "Hmax", 10, "Hmin", 9);
%model = generateMesh(femodel);

% Show the mesh
%pdemesh(model)

%%

Nodes    = model.Geometry.Mesh.Nodes./1000;% Scale from [mm] to [m]. The mesh should be already with [m] units!
Elements = model.Geometry.Mesh.Elements;

%%
% Build the LVPBody

Primitives = {PCStretchCompressionPrimitive(0.3), ...
              PCTwistPrimitive(0.3), ...
              PCShearPrimitive(0.3), ...
              PCBendingPrimitive(0.3)};

Primitives = {PCStretchCompressionPrimitive(0.3), ...
              PCTwistShearPrimitive(0.3), ...
              PCBendingPrimitive(0.3)};

B = LVPBody(Nodes, Elements, Primitives, 1);

%% Do a plot
close all;
figure; hold on; grid on; view(3)
light("Position", [-0.1, -0.1, 0.1])
lighting gouraud
%q_test = [0; 0; 0; 0; 0; pi];

q_test = [-1; 10; -0.1; 0.1; 2; 1];

B.plot(q_test, "LineStyle", "-", "FaceAlpha", 1);

axis equal

%% Test code generation
codegen -o test/test_lvpbody/lvp_codegen_mex lvp_codegen -args {Nodes, Elements}





