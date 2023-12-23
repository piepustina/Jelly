%% Create an instance of the class
clear; clc;
DHTable         = [0 1 0 0; 
                  -1 0 1 0];
m               = 1;
I               = eye(3);
pcom            = [1; 0; -1];
g               = [-9.81; 0; 0];
n               = size(DHTable, 1);
N_B             = n;
B1 = cell(N_B, 1);
J1 = cell(N_B, 1);
for i = 1:N_B
    Parameters  = [m; pcom; I(1, 1); I(2, 2); I(3, 3); I(1, 2); I(1, 3); I(2, 3)];
    B1{i}       = RigidBody(Parameters);
    J1{i}       = RotationalJoint(DHTable(i, :)');
end

%Build an equivalent body for comparison
b2 = PCCCylindricalDelta_Body(Parameters);
j2 = PCCCylindricalDelta_Joint(Parameters(1));
B2 = cell(N_B, 1);
J2 = cell(N_B, 1);
for i = 1:N_B
    Parameters  = [DHTable(i, 1:3)'; m; pcom; I(1, 1); I(2, 2); I(3, 3); I(1, 2); I(1, 3); I(2, 3)];
    B2{i} = Rotational_Body(Parameters);
    J2{i} = Rotational_Joint(Parameters);
end
%Test configurations
q_test      = 1;
dq_test     = 1;
ddq_test    = 1;


%% Test dynamic methods

r1 = BodyTree(J1, B1);
r1.g = g;
r2 = BodyTreeOld(J2, B2, 1e-4*ones(n*N_B, 1), 1e-4*ones(n*N_B, 1), 1e-4*ones(n*N_B, 1), 1e-4*ones(n*N_B, 1));
r2.g = g;

%% Numerical test
disp("Mass matrix...")
r1.MassMatrix(repmat(q_test, N_B, 1))
r2.MassMatrix(repmat(q_test, N_B, 1), 'double')


disp("Gravity force...")
r1.GravityForce(repmat(q_test, N_B, 1))
r2.GravityForce(repmat(q_test, N_B, 1), 'double')

%% Symbolic test
q  = sym('q', [n, 1], 'real');
dq = sym('dq', [n, 1], 'real');
r1.MassMatrix(q)
r2.MassMatrix(q, 'sym')

r1.GravityForce(q)
r2.GravityForce(q, 'sym')


%%

rStruct = TreeStructConverter.ObjectToStruct(r1);

r3      = TreeStructConverter.StructToObject(rStruct);



