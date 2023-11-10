%Run this script to compute the model of a PCC planar slender soft body and
%of its joint. This script uses the cosine function representation of the
%sinc function.
%Initialize the toolbox
initToolbox;

tic;
%Compute all the model of a PCC planar without elongation slender body and of its joint
%Define the configuration and its time derivatives
syms theta dtheta ddtheta real

%Define the parameters of the body (mass density and rest length)
%In practice, we assume to know the mass and the rest length
syms rho L_0 E Poi Eta R real
%E = Young's modulus
%Poi = Poisson ration
%Eta = Material damping
%R   = cross sectional radius
%Rest length
L = L_0;

%Mass densiy along the curvilinar abscissa
rho_L = rho*pi*R^2;

%Mass of the segment
m = rho_L*L;

%
%Inertia of each body in local frame (possibly as a function of L)
I_body = rho*diag([pi*R^4/4, pi*R^4/4, pi*R^4/2]);


q   = [theta];
dq  = [dtheta];
ddq = [ddtheta];

params = [L_0; R; rho; E; Poi; Eta];

%Define the variables of integration
syms s real

syms x real
syms SINC(x) COSINC(x)

SINC(x)   = 1/(gamma(1-x/pi)*gamma(1+x/pi));
k = 3;
SINC(x)   = prod(arrayfun(@(i) cos(x/(2^i)), 1:k));
COSINC(x) = SINC(x/2)*sin(x/2);


%Use the simple definition of the sinc and cosinc for the strains
T_s_strain = [cos(s/L*theta), 0, -sin(s/L*theta), -L*(1-cos(s/L*theta))/theta;
       0             , 1,               0, 0;
       sin(s/L*theta), 0,  cos(s/L*theta), L*sin(s/L*theta)/theta    ;
              0,       0, 0,                      1];
%Use the gamma definition of the sinc and cosinc for the kinematics
T_s = [cos(s/L*theta), 0, -sin(s/L*theta), -COSINC(s/L*theta)*s;
       0             , 1,               0, 0;
       sin(s/L*theta), 0,  cos(s/L*theta), SINC(s/L*theta)*s; 
             0       , 0,               0,                  1];
T_s = simplify(T_s);
T = subs(T_s, s, L);

T_inv = simplify(invTransformation(T));            
disp("The transformation matrix is:");
pretty(T)

%% Compute the strain field
xi_hat = simplify(T_s_strain\diff(T_s_strain, s));

%Extract the strain as a 6D vector in angular and linear part
xi = [skew_inv(xi_hat(1:3, 1:3)); xi_hat(1:3, 4)];

%% Compute the kinematic vectors
%Take the rotation matrix and the position vector
R_i_1_i = T(1:3, 1:3);
t_i_1_i = T(1:3, 4);

prev_p_i = T_s(1:3, 4);

%Rotate the vectors in {S_i}
disp("Computing p_i_s...");
p_i    = T_inv*[prev_p_i;1];
p_i    = simplify(p_i(1:3), 'Seconds', 5);
disp("Computing dp_i_s...");
dp_i   = matrix_derivative(p_i, q, dq);
dp_i   = simplify(dp_i, 'Seconds', 5);
disp("Computing ddp_i_s...");
ddp_i  = matrix_derivative(dp_i, [q; dq], [dq; ddq]);
ddp_i   = simplify(ddp_i, 'Seconds', 5);


p_comi_s   = (1/m)*p_i*rho_L;
dp_comi_s  = (1/m)*dp_i*rho_L;
ddp_comi_s = (1/m)*ddp_i*rho_L;

%% Compute all the required terms

className = "PCCPlanar_Body_no_elongation_numerical";
%Compute all the required terms
destination_folder_j = "./classes/@PCCPlanar_Joint_no_elongation_numerical";
destination_folder_b = "./classes/@PCCPlanar_Body_no_elongation_numerical";

disp("Computing the data of the model...");
hold_integrals = true;
IntRelTol = 1e-3;
IntAbsTol = 1e-4;
compute_model_data_numerically;

%% Export the functions

export_numerical_functions;