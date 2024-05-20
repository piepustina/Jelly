% Load a test mesh
clc; clear; close all;
%femodel = femodel(Geometry="./test/meshes/Cylinder_coarse.stl");
femodel = femodel(Geometry="./test/meshes/Diamond.stl");


%model = generateMesh(femodel, "Hmax", 60, "Hmin", 60, "Hgrad", 2);
model = generateMesh(femodel, "Hmax", 10, "Hmin", 9);
%model = generateMesh(femodel);

% Show the mesh
%pdemesh(model)

%%

Nodes    = model.Geometry.Mesh.Nodes./1000;% Scale from [mm] to [m]. The mesh should be already with [m] units!
Elements = model.Geometry.Mesh.Elements;

%%
% Build the LVPBody

L0 = max(Nodes(3, :));
Primitives = {PCStretchCompressionPrimitive(L0), ...
              PCTwistShearPrimitive(L0), ...
              PCBendingPrimitive(L0)};

B = LVPBody(Nodes, Elements, Primitives, 1);

B1 = LVPBody(Nodes, Elements, Primitives, 1);

%% Do a plot

r = BodyTree({FixedJoint(); RotationalJoint(zeros(4, 1))}, {B; B1});

close all;
figure; hold on; grid on; view(3)
light("Position", [-0.1, -0.1, 0.1])
lighting gouraud

q_test = 0*[-0.04; 1; -0.01; 0.01; 0.2; 0.1; 0; 0; 1; -0.01; 0.01; 0.2; 0.1];

r.plot(q_test, "LineStyle", "-", "FaceAlpha", 1);

xlabel("$x$", "Interpreter", "latex", "FontSize", 14)
ylabel("$y$", "Interpreter", "latex", "FontSize", 14)
zlabel("$z$", "Interpreter", "latex", "FontSize", 14)

axis equal
view(3)

%% Test code generation
codegen -o test/test_lvpbody/lvp_codegen_mex lvp_codegen -args {Nodes, Elements}

%% Test code generation with cuda
codegen -config coder.gpuConfig('mex') -o test/test_lvpbody/lvp_codegen_mex_cuda lvp_codegen -args {Nodes, Elements}
