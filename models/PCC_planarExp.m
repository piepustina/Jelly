%Run this script to compute the model of a PCC planar slender soft body and
%of its joint. 
%Initialize the toolbox
initToolbox;

%Compute all the model of a PCC planar slender body and of its joint.

%Define the configuration and its time derivatives
syms theta deltaL dtheta ddeltaL ddtheta dddeltaL real

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


q   = [theta; deltaL];
dq  = [dtheta; ddeltaL];
ddq = [ddtheta; dddeltaL];

params = [L_0; R; rho; E; Poi; Eta];

%Define the variables of integration
syms s real

%%
%Extract the strain as a 6D vector in angular and linear part
xi       = [0; 0; q(1)/L_0; (L_0+q(2))/L_0; 0; 0];
xi_hat   = [skew(xi(1:3)), xi(4:6); zeros(1, 4)];

simplify(expm(xi_hat*s))

%%
%Compute the transformation matrix
Omega = Omega_1 + Omega_2 + Omega_3 + Omega_4;

%Compute the matrix exponential
Theta = simplify(sqrt(skew_inv(Omega(1:3, 1:3))'*skew_inv(Omega(1:3, 1:3))), 'Seconds', 10);
disp("Computing T_s");

T_s   = (eye(4) + Omega + ...
          1/(Theta^2)*(1-cos(Theta))*mpower(Omega, 2) +...
          1/(Theta^3)*(Theta-sin(Theta))*mpower(Omega, 3))';

%%

T = subs(T_s, s, L);
            
T_inv = simplify(invTransformation(T));            
disp("The transformation matrix is:");
pretty(T)

%Take the rotation matrix and the position vector
R_i_1_i = T(1:3, 1:3);
t_i_1_i = T(1:3, 4);

%
prev_p_i = T_s(1:3, 4);

%Compute the CoM in {S_{i-1}}
prev_p_comi = (1/m)*int(prev_p_i*rho_L, s, [0, L]);


prev_p_comi = simplify(prev_p_comi);
prev_r_i    = simplify(prev_p_i - prev_p_comi);


%Rotate the vectors in {S_i}
p_comi = T_inv*[prev_p_comi;1];
p_i    = T_inv*[prev_p_i;1];

p_comi = simplify(p_comi(1:3));
p_i    = simplify(p_i(1:3));

%Compute the relative position vector and its derivatives
r_i    = p_i - p_comi;
dr_i   = matrix_derivative(r_i, q, dq);
ddr_i  = matrix_derivative(dr_i, [q;dq], [dq;ddq]);

assert(~all(simplify(int(r_i*rho_L, s, [0, L]))))

%%

%Compute all the required terms
%Set to true to leave the integrals in implicit form
hold_integrals = false;

destination_folder_j = "./classes/@PCCPlanar_Joint";
destination_folder_b = "./classes/@PCCPlanar_Body";

compute_model_data;

disp("Export the functions...");
export_functions;