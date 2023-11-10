%Compute the model of a piecewise affine body with cylindrical shape.
%All the inertial parameters, e.g., mass density, Young's modulus, etc...,
%are assumed constant.
%Make sure that the toolbox has been initialized

clear, clc;

%Define the configuration and its time derivatives (no shear)
syms kappa_X kappa_Y tau_Z gamma_X gamma_Y deltaL ...
     dkappa_X dkappa_Y dtau_Z dgamma_X dgamma_Y ddeltaL...
     ddkappa_X ddkappa_Y ddtau_Z ddgamma_X ddgamma_Y dddeltaL real

% L_0 = Rest length
% rho = mass density of the body
% R   = cross sectional radius
% E = Young's modulus
% Poi = Poisson ration
% Eta = Material damping

syms L_0 R rho E_mod Poi Eta real positive

E = E_mod;

%Reference length used to compute the integral
L = L_0;
%Mass density along the curvilinar abscissa
rho_L = rho*pi*R^2;

%Mass of the segment
m = rho_L*L_0;

%
q   = [kappa_X; kappa_Y; tau_Z; gamma_X; gamma_Y; deltaL];
dq   = [dkappa_X; dkappa_Y; dtau_Z; dgamma_X; dgamma_Y; ddeltaL];
ddq   = [ddkappa_X; ddkappa_Y; ddtau_Z; ddgamma_X; ddgamma_Y; dddeltaL];

params = [L_0; R; rho; E; Poi; Eta];

syms s real positive

%% Define the strain field
%xi_hat = simplify(T_s\diff(T_s, s));

%Extract the strain as a 6D vector in angular and linear part
xi       = [-q(2)/L; q(1)/L; q(3)/L; q(4)/L; q(5)/L; (L+q(6))/L];
xi_hat   = [skew(xi(1:3)), xi(4:6); zeros(1, 4)];
xi_s     = s*xi;
xi_hat_s = s*xi_hat;

%Compute the matrix exponential, i.e., exp(xi_hat_s)
Theta = sqrt(xi_s(1:3)'*xi_s(1:3));
T_s   = simplify(eye(4) + xi_hat_s + ...
          1/(Theta^2)*(1-cos(Theta))*mpower(xi_hat_s, 2) +...
          1/(Theta^3)*(Theta-sin(Theta))*mpower(xi_hat_s, 3));
% T_s = expm(xi_hat_s);

T  = subs(T_s, s, L);

T = simplify(T);
%%
T_inv = simplify(invTransformation(T));            
disp("The transformation matrix is:");
pretty(T)

%Take the rotation matrix and the position vector
R_i_1_i = T(1:3, 1:3);
t_i_1_i = T(1:3, 4);

prev_p_i = T_s(1:3, 4);

%Compute the CoM in {S_{i-1}}
% prev_p_comi = (1/m)*int(prev_p_i*rho_L, s, [0, L]);
% prev_p_comi = simplify(prev_p_comi);
% prev_r_i    = simplify(prev_p_i - prev_p_comi);

%Rotate the vectors in {S_i}
disp("Computing p_i_s...");
%p_i    = T_inv*[prev_p_i;1];
p_i  = R_i_1_i'*( prev_p_i - subs(prev_p_i, s, L));
p_i    = simplify(p_i(1:3), 'Seconds', 10);
disp("Computing dp_i_s...");
dprev_p_i = matrix_derivative(prev_p_i, q, dq);
dR_i_1_i = matrix_derivative(R_i_1_i, q, dq);
dp_i = dR_i_1_i'*( prev_p_i - subs(prev_p_i, s, L)) +...
        R_i_1_i'*(dprev_p_i - subs(dprev_p_i, s, L));
dp_i   = simplify(dp_i, 'Seconds', 10);
disp("Computing ddp_i_s...");
ddprev_p_i = matrix_derivative(dprev_p_i, [q; dq], [dq; ddq]);
ddR_i_1_i  = matrix_derivative(dR_i_1_i,  [q; dq], [dq; ddq]);
ddp_i      = ddR_i_1_i'*(prev_p_i - subs(prev_p_i, s, L)) +...
      +2*dR_i_1_i'*( dprev_p_i - subs(dprev_p_i, s, L)) +...
         R_i_1_i'*(ddprev_p_i - subs(ddprev_p_i, s, L));
ddp_i   = simplify(ddp_i, 'Seconds', 10);

p_comi_s   = (1/m)*p_i*rho_L;
dp_comi_s  = (1/m)*dp_i*rho_L;
ddp_comi_s = (1/m)*ddp_i*rho_L;

% p_comi = int(p_comi_s, s, [0, L]);
% p_comi = simplify(p_comi, 'Seconds', 10);
% 
% % %Compute the relative position vector and its derivatives
% r_i    = p_i - p_comi;
% dr_i   = matrix_derivative(r_i, q, dq);
% ddr_i  = matrix_derivative(dr_i, [q;dq], [dq;ddq]);
% 
% r_i = simplify(r_i, 'Seconds', 10);
% dr_i = simplify(dr_i, 'Seconds', 10);
% ddr_i = simplify(ddr_i, 'Seconds', 10);

%assert(~all(simplify(int(r_i*rho_L, s, [0, L]))))

%% Compute all the required terms

className = "PCSCylindrical_Body";
destination_folder_b = "./classes/@PCSCylindrical_Body";
destination_folder_j = "./classes/@PCSCylindrical_Joint";

disp("Computing the data of the model...");
hold_integrals = true;
IntRelTol = 1e-3;
IntAbsTol = 1e-4;
compute_model_data_numerically;
%compute_model_data;

%% Export the functions
export_numerical_functions;
%export_functions;
