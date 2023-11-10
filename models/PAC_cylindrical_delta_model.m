%Compute the model of a piecewise affine body with cylindrical shape.
%All the inertial parameters, e.g., mass density, Young's modulus, etc...,
%are assumed constant.
%Make sure that the toolbox has been initialized

clear, clc;

%Define the configuration and its time derivatives
syms theta_X_1 theta_X_2 theta_Y_1 theta_Y_2 deltaL ...
     dtheta_X_1 dtheta_X_2 dtheta_Y_1 dtheta_Y_2 ddeltaL...
     ddtheta_X_1 ddtheta_X_2 ddtheta_Y_1 ddtheta_Y_2 dddeltaL real

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
q   = [theta_X_1; theta_X_2; theta_Y_1; theta_Y_2; deltaL];
dq  = [dtheta_X_1; dtheta_X_2; dtheta_Y_1; dtheta_Y_2; ddeltaL];
ddq = [ddtheta_X_1; ddtheta_X_2; ddtheta_Y_1; ddtheta_Y_2; dddeltaL];

params = [L_0; R; rho; E; Poi; Eta];

syms s real positive

%% Define the strain field
%xi_hat = simplify(T_s\diff(T_s, s));

%Extract the strain as a 6D vector in angular and linear part
xi       = [-q(3)/L - q(4)*s/L^2; q(1)/L + q(2)*s/L^2; 0; 0; 0; (L+q(5))/L];
xi_hat   = [skew(xi(1:3)), xi(4:6); zeros(1, 4)];

xi_hat   = [skew(xi(1:3)), xi(4:6); zeros(1, 4)]';

%Compute the Magnus approximation of the strain field
syms t1 t2 t3 t4 real
xi_hat_t1 = subs(xi_hat, s, t1);
xi_hat_t2 = subs(xi_hat, s, t2);
xi_hat_t3 = subs(xi_hat, s, t3);
xi_hat_t4 = subs(xi_hat, s, t4);
Omega_1 = int(xi_hat_t1, t1, [0, s]);
Omega_2 = 1/2*int(int(commutator(xi_hat_t1, xi_hat_t2), t2, [0, t1]), t1, [0, s]);
integrand_3 = commutator(xi_hat_t1, commutator(xi_hat_t2, xi_hat_t3)) + ...
              commutator(xi_hat_t3, commutator(xi_hat_t2, xi_hat_t1));
Omega_3 = 1/6*int(int(int(integrand_3, t3, [0, t2]), t2, [0, t1]), t1, [0, s]);
integrand_4 = commutator(commutator(commutator(xi_hat_t1, xi_hat_t2), xi_hat_t3), xi_hat_t4) +...
              commutator(xi_hat_t1, commutator(commutator(xi_hat_t2, xi_hat_t3), xi_hat_t4)) + ...
              commutator(xi_hat_t1, commutator(xi_hat_t2, commutator(xi_hat_t3, xi_hat_t4))) + ...
              commutator(xi_hat_t2, commutator(xi_hat_t3, commutator(xi_hat_t4, xi_hat_t1)));
Omega_4 = 1/12*int(int(int(int(integrand_4, t4, [0, t3]), t3, [0, t2]), t2, [0, t1]), t1, [0, s]);
Omega = Omega_1 + Omega_2 + Omega_3 + Omega_4;

%Compute the matrix exponential
Theta = simplify(sqrt(skew_inv(Omega(1:3, 1:3))'*skew_inv(Omega(1:3, 1:3))), 'Seconds', 10);
disp("Computing T_s");
% T_s   = simplify(eye(4) + Omega + ...
%           1/(Theta^2)*(1-cos(Theta))*mpower(Omega, 2) +...
%           1/(Theta^3)*(Theta-sin(Theta))*mpower(Omega, 3), 'Seconds', 10);

T_s   = (eye(4) + Omega + ...
          1/(Theta^2)*(1-cos(Theta))*mpower(Omega, 2) +...
          1/(Theta^3)*(Theta-sin(Theta))*mpower(Omega, 3))';

% T_s1_ang = eye(3)+sin(Theta)/Theta*Omega(1:3, 1:3)+(1-cos(Theta))/(Theta^2)*mpower(Omega(1:3, 1:3), 2);
% T_s1_lin = eye(3)+(1-cos(Theta))/Theta^2*Omega(1:3, 1:3)+(Theta-sin(Theta))/(Theta^3)*mpower(Omega(1:3, 1:3), 2);
% T_s1_lin = T_s1_lin*Omega(1:3, 4);


T  = subs(T_s, s, L);

disp("Simplyfing T");
T = simplify(T, 'Seconds', 10);

disp("Inverting T");
T_inv = simplify(invTransformation(T), 'Seconds', 10);            
disp("The transformation matrix is:");
pretty(T)

%%
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

% p_comi = int(p_comi_s, s, [0, L]);
% p_comi = simplify(p_comi);

% %Compute the relative position vector and its derivatives
% r_i    = p_i - p_comi;
% dr_i   = matrix_derivative(r_i, q, dq);
% ddr_i  = matrix_derivative(dr_i, [q;dq], [dq;ddq]);
% 
% r_i = simplify(r_i);
% dr_i = simplify(dr_i);
% ddr_i = simplify(ddr_i);

%assert(~all(simplify(int(r_i*rho_L, s, [0, L]))))

%% Compute all the required terms

className = "PACCylindrical_Delta_Body";
destination_folder_b = "./classes/@PACCylindrical_Delta_Body";
destination_folder_j = "./classes/@PACCylindrical_Delta_Joint";

disp("Computing the data of the model...");
hold_integrals = true;
IntRelTol = 1e-3;
IntAbsTol = 1e-4;
compute_model_data_numerically;
%compute_model_data;

%% Export the functions

export_numerical_functions;
%export_functions;
