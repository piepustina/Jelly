%%
clear; clc;
initToolbox;

%%
n_B = 2;
L_0 = sym('L0', [n_B, 1], {'real', 'positive'});
R   = sym('R', [n_B, 1], {'real', 'positive'});
rho = sym('rho', [n_B, 1], {'real', 'positive'});
E   = sym('E', [n_B, 1], {'real', 'positive'});
Poi = sym('Poi', [n_B, 1], {'real', 'positive'});
Eta = sym('Eta', [n_B, 1], {'real', 'positive'});

%Define tread-like actuators
D  = sym('D', [1, 1], {'real', 'positive'});
Actuators = {};
Actuators{1} = ConstantActuator(D*[1; 0; 0], [0 L_0(1)+L_0(2)], [1 2]);
Actuators{2} = ConstantActuator(D*[cos(2*pi/3); sin(2*pi/3); 0], [0 L_0(1)+L_0(2)], [1 2]);
Actuators{3} = ConstantActuator(D*[cos(4*pi/3); sin(4*pi/3); 0], [0 L_0(1)+L_0(2)], [1 2]);

% 
% Actuators{1} = ObliqueActuator(D*[0 0; 1 0; 0 0], [0 L_0(1)+L_0(2)], [1 2]);
% Actuators{2} = ObliqueActuator(D*[sin(2*pi/3) 0; cos(2*pi/3) 0; 0 0], [0 L_0(1)+L_0(2)], [1 2]);
% Actuators{3} = ObliqueActuator(D*[sin(4*pi/3) 0; cos(4*pi/3) 0; 0 0], [0 L_0(1)+L_0(2)], [1 2]);


%Define the joints of the robot
Joints = cell(n_B, 1);
Bodies = cell(n_B, 1);

for i = 1:n_B
    Parameters = [L_0(i), R(i), rho(i), E(i), Poi(i), Eta(i)];
    Joints{i} = PCCCylindrical_Joint(Parameters(1));
    Bodies{i} = PCCCylindrical_Body(Parameters);
end

%Construct a BodyTree without computing the actuation matrix
bt_3D = BodyTree(Joints, Bodies, Actuators, false);

%%
%Compute the actuation matrix symbolically
clc
[q, Aq, y] = bt_3D.ComputeActuationMatrix();
Aq = simplify(Aq);

%Sanity check that Aq is the Jacobian of y
simplify(jacobian(y, q)' - Aq) 

%%
%Define the change of coordinates
theta_ = [y; q(1:3)];

theta_ = [y;
          q(3) + D*(q(1)*cos(q(2))); 
          q(3) - D/2*q(1)*cos(q(2)) + sqrt(3)/2*D*q(1)*sin(q(2)); 
          q(3) - D/2*q(1)*cos(q(2)) - sqrt(3)/2*D*q(1)*sin(q(2))];

%Check that theta qualifies as a change of coordinates
jTheta = jacobian(theta_, q);
rank(jTheta)

%Compute the new actuation matrix
simplify(jTheta'\Aq)

%% Check that this is correct
matlabFunction(theta_, 'Vars', {q, D}, 'File', 'trash/theta_f');

addpath('trash')
