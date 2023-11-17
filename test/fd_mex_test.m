function dx = fd_mex_test(x)
%FD_MEX_TEST Test code generation on the forward dynamics

DHTable         = [0 1 0 0; 
                  -1 0 1 0];
m               = 1;
I               = eye(3);
pcom            = [1; 0; -1];
g               = [-9.81; 0; 0];
n               = size(DHTable, 1);
N_B             = n;
B1 = cell(BodyTree.MaxBodiesNumber, 1);
J1 = cell(BodyTree.MaxBodiesNumber, 1);
for i = 1:BodyTree.MaxBodiesNumber
    if i <= N_B
        Parameters  = [m; pcom; I(1, 1); I(2, 2); I(3, 3); I(1, 2); I(1, 3); I(2, 3)];
        B1{i}       = RigidBody(Parameters);
        J1{i}       = RotationalJoint(DHTable(i, :)');
    else
        B1{i} = 0;
        J1{i} = 0;
    end
end

%Create the body tree
r1 = BodyTree(J1, B1);
r1.g = g;

%Compute the forward dynamics
dx = r1.StateSpaceForwardDynamics(0, x, zeros(r1.n, 1));

end

