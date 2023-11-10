%% Perform the computations of the integrals numerically
%--------------------------------------------------------------------------
% Compute relative velocities and accelerations of the planar PCC joint----
%--------------------------------------------------------------------------
disp("Computing kinematic data...");
disp("Relative linear velocity...");
%Relative linear velocity
v_i_1_i = matrix_derivative(t_i_1_i, q, dq);
v_i_1_i = simplify(v_i_1_i, 'Seconds', 10);
%Relative angular velocity
disp("Relative angular velocity...");
S_omega_i_1_i = matrix_derivative(R_i_1_i, q, dq)*R_i_1_i';
omega_i_1_i   = [S_omega_i_1_i(3,2);S_omega_i_1_i(1,3);S_omega_i_1_i(2,1)];
omega_i_1_i = simplify(omega_i_1_i, 'Seconds', 10);
%Relative linear acceleration
disp("Relative linear acceleration...");
dv_i_1_i = matrix_derivative(v_i_1_i, [q;dq], [dq;ddq]);
%Relative angular acceleration
disp("Relative angular acceleration...");
domega_i_1_i = matrix_derivative(omega_i_1_i, [q;dq], [dq;ddq]);
%Relative CoM velocity
disp("Relative CoM velocity...");
%v_i_comi = matrix_derivative(p_comi, q, dq);
v_i_comi = str2func(sprintf("@(q, dq, params) myintegral(@(s) %s.dp_comi_s(q, dq, params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
            className, IntAbsTol, IntRelTol));
%Relative CoM acceleration
disp("Relative CoM acceleration...");
%dv_i_comi = matrix_derivative(v_i_comi, [q; dq], [dq; ddq]);
dv_i_comi = str2func(sprintf("@(q, dq, ddq, params) myintegral(@(s) %s.ddp_comi_s(q, dq, ddq, params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
            className, IntAbsTol, IntRelTol));
%Partial linear velocity
disp("Partial linear velocity...");
v_par_i   = R_i_1_i'*jacobian(v_i_1_i, dq);
%Partial angular velocity
disp("Partial angular velocity...");
omega_par_i   = R_i_1_i'*jacobian(omega_i_1_i, dq);
%%
disp("Computing remaining terms...");
%Strain field
xi_i = xi;
%Strain field in the reference configuration
xi_i_star = subs(xi_i, q, zeros(size(q)));
%Partial strain
xi_par_i  = jacobian(xi, q);
%Jacobian of the relative position vector computed differently
%J_dr_i = jacobian(r_i, q);
J_p_i      = jacobian(p_i, q);
%J_p_i      = simplify(J_p_i, 'Seconds', 10);
%J_p_comi_s = jacobian(p_comi_s, q);
J_p_comi_s = (1/m)*J_p_i*rho_L;
J_dr_i = J_p_i - int(J_p_comi_s, s, [0, L], 'Hold', hold_integrals);

%% Approximate numerically all the required integrals
%Export p_i and its time derivatives
disp("Exporting p_i and derivatives...");
matlabFunction(p_i,  'Vars', {q, params, s}, 'File', destination_folder_b + "/p_i_s");
matlabFunction(dp_i, 'Vars', {q, dq, params, s}, 'File', destination_folder_b + "/dp_i_s");
matlabFunction(ddp_i,'Vars', {q, dq, ddq, params, s}, 'File', destination_folder_b + "/ddp_i_s");

% matlabFunction(p_comi_s, 'Vars', {q, params, s}, 'File',  destination_folder_b + "/p_comi_s", 'Optimize', true);
% matlabFunction(dp_comi_s, 'Vars', {q, dq, params, s}, 'File',  destination_folder_b + "/dp_comi_s", 'Optimize', true);
% matlabFunction(ddp_comi_s, 'Vars', {q, dq, ddq, params, s}, 'File',  destination_folder_b + "/ddp_comi_s", 'Optimize', true);

p_comi_s_a = str2func(sprintf("@(q, params, s) (1/%s.m(params))*(%s.p_i_s(q, params, s))*(%s.rho_L_s(params, s))", className, className, className));
dp_comi_s_a = str2func(sprintf("@(q, dq, params, s) (1/%s.m(params))*(%s.dp_i_s(q, dq, params, s))*(%s.rho_L_s(params, s))", className, className, className));
ddp_comi_s_a = str2func(sprintf("@(q, dq, ddq, params, s) (1/%s.m(params))*(%s.ddp_i_s(q, dq, ddq, params, s))*(%s.rho_L_s(params, s))", className, className, className));

p_comi = str2func(sprintf("@(q, params) myintegral(@(s) %s.p_comi_s(q, params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", className, IntAbsTol, IntRelTol));
r_i_s   = str2func(sprintf("@(q, params, s) %s.p_i_s(q, params, s)-myintegral(@(s) %s.p_comi_s(q, params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
            className, className, IntAbsTol, IntRelTol));
dr_i_s  = str2func(sprintf("@(q, dq, params, s) %s.dp_i_s(q, dq, params, s)-myintegral(@(s) %s.dp_comi_s(q, dq, params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
            className, className, IntAbsTol, IntRelTol));
ddr_i_s = str2func(sprintf("@(q, dq, ddq, params, s) %s.ddp_i_s(q, dq, ddq, params, s)-myintegral(@(s) %s.ddp_comi_s(q, dq, ddq, params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
            className, className, IntAbsTol, IntRelTol));
matlabFunction(rho_L, 'Vars', {params, s}, "File", destination_folder_b + "/rho_L_s");
matlabFunction(J_p_i, 'Vars', {q, params, s}, 'File', destination_folder_b + "/J_p_i");
matlabFunction(J_p_comi_s, 'Vars', {q, params, s}, 'File', destination_folder_b + "/J_p_comi_s");
J_dr_i_s = str2func(sprintf("@(q, params, s) %s.J_p_i(q, params, s) - myintegral(@(s) %s.J_p_comi_s(q, params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
            className, className, IntAbsTol, IntRelTol));
%J_r_i_X_dr_i_s = matlabFunction(skew(r_i)*J_dr_i*rho_L, 'Vars', {q, params, s}, "File", destination_folder_b + "/J_r_i_X_dr_i_s");
J_r_i_X_dr_i_s = str2func(sprintf("@(q, params, s) skew(%s.r_i_s(q, params, s))* %s.J_dr_i_s(q, params, s)*%s.rho_L_s(params, s)", ...
                    className, className, className));

%%

%--------------------------------------------------------------------------
% Compute all the needed integrals ----------------------------------------
%--------------------------------------------------------------------------
disp("Computing integrals...");
disp("Computing I");
I = str2func(sprintf("@(q, params) myintegral(@(s) skew(%s.r_i_s(q, params, s))'*skew(%s.r_i_s(q, params, s))*%s.rho_L_s(params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
            className, className, className, IntAbsTol, IntRelTol));

%--------------------------------------------------------------------------
disp("Computing dI");
dI = zeros(3, 3, 'sym');

%--------------------------------------------------------------------------
disp("Computing J");
J = str2func(sprintf("@(q, dq, params) myintegral(@(s) skew(%s.dr_i_s(q, dq, params, s))'*skew(%s.r_i_s(q, params, s))*%s.rho_L_s(params, s) + skew(%s.r_i_s(q, params, s))'*skew(%s.dr_i_s(q, dq, params, s))*%s.rho_L_s(params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
    className, className, className, className, className, className, IntAbsTol, IntRelTol));

%--------------------------------------------------------------------------
disp("Computing int_dr_i");
int_dr_i = zeros(3, 1, 'sym');

%--------------------------------------------------------------------------
disp("Computing int_ddr_i");
int_ddr_i = zeros(3, 1, 'sym');

%--------------------------------------------------------------------------
disp("Computing int_r_i_X_dr_i");

int_r_i_X_dr_i = str2func(sprintf("@(q, dq, params) myintegral(@(s) cross(%s.r_i_s(q, params, s), %s.dr_i_s(q, dq, params, s))*%s.rho_L_s(params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
    className, className, className, IntAbsTol, IntRelTol));

%--------------------------------------------------------------------------
disp("Computing int_r_i_X_ddr_i");
int_r_i_X_ddr_i = str2func(sprintf("@(q, dq, ddq, params) myintegral(@(s) cross(%s.r_i_s(q, params, s), %s.ddr_i_s(q, dq, ddq, params, s))*%s.rho_L_s(params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
    className, className, className, IntAbsTol, IntRelTol));

%--------------------------------------------------------------------------
disp("Computing int_dr_i_X_pv_r");
int_dr_i_X_pv_r = str2func(sprintf("@(q, dq, params) myintegral(@(s) %s.J_dr_i_s(q, params, s)'*skew(%s.dr_i_s(q, dq, params, s))'*%s.rho_L_s(params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
    className, className, className, IntAbsTol, IntRelTol));

%--------------------------------------------------------------------------
disp("Computing int_pv_r_O_dd_r");

int_pv_r_O_dd_r = str2func(sprintf("@(q, dq, ddq, params) myintegral(@(s) %s.J_dr_i_s(q, params, s)'*%s.ddr_i_s(q, dq, ddq, params, s)*%s.rho_L_s(params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
    className, className, className, IntAbsTol, IntRelTol));

%--------------------------------------------------------------------------
disp("Computing int_dr_i_O_dr_i");
int_dr_i_O_dr_i = zeros([1,1], 'sym');

%%
%--------------------------------------------------------------------------
% Compute all the needed gradients ----------------------------------------
%--------------------------------------------------------------------------
disp("Computing gradients...");
grad_int_dr_i = zeros(size(J_dr_i'), 'sym');

grad_int_r_i_X_dr_i = str2func(sprintf("@(q, params) myintegral(@(s) %s.J_r_i_X_dr_i_s(q, params, s)', 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
                              className, IntAbsTol, IntRelTol));

%grad_J_s_sym = (sympagemtimes(pageskew(J_dr_i, 1), skew(r_i))+sympagemtimes(skew(r_i)', pageskew(J_dr_i, 0)))*rho_L
%grad_J_s = matlabFunction(grad_J_s_sym, 'Vars', {q, params, s}, "File", destination_folder_b + "/grad_J_s");
grad_J_s = str2func(sprintf("@(q, params, s) (sympagemtimes(pageskew(%s.J_dr_i_s(q, params, s), 1), skew(%s.r_i_s(q, params, s)))+sympagemtimes(skew(%s.r_i_s(q, params, s))', pageskew(%s.J_dr_i_s(q, params, s), 0)))*%s.rho_L_s(params, s)", ...
                            className, className, className, className, className));

grad_J   = str2func(sprintf("@(q, params) myintegral(@(s) %s.grad_J_s(q, params, s), 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
            className, IntAbsTol, IntRelTol));
%grad_v_com_i = jacobian(v_i_comi, dq)' = int(J_p_comi_s', s, [0, L]);
grad_v_com_i = str2func(sprintf("@(q, params) myintegral(@(s) %s.J_p_comi_s(q, params, s)', 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", ...
            className, IntAbsTol, IntRelTol));

%%
%--------------------------------------------------------------------------
% Compute the generalized internal forces ---------------------------------
%--------------------------------------------------------------------------
disp("Computing internal forces...");

% Elastic force
%Shear modulus (material is assumed isotropic)
G = E/(2*(1+Poi));
%Inertia of each body in local frame for the stiffness
I_body = diag([pi*R^4/4, pi*R^4/4, pi*R^4/2]);
%Compute the linear and angular stiffness
K_l = pi*R^2*diag([G, G, E]);
K_a = diag([E*I_body(1, 1), E*I_body(2, 2), G*I_body(3, 3)]);
%Body stiffness
K_b = blkdiag(K_a, K_l);
%Body elastic force
dF_elastica = K_b*(xi_i - xi_i_star);

% Damping force
D_l = K_l*Eta;
D_a = K_a*Eta;
%Body damping matrix
D_b = blkdiag(D_a, D_l);

%Body damping force
dF_damping = D_b*matrix_derivative(xi_i, q, dq);

%Computation of the generalized forces by projecting the forces into the
%strain space
K_s   = xi_par_i'*dF_elastica;
matlabFunction(K_s, 'Vars', {q, params, s}, "File", destination_folder_b + "/K_s");
K_q = str2func(sprintf("@(q, params) myintegral(@(s) %s.K_s(q, params, s)  , 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", className, IntAbsTol, IntRelTol));
D_s   = xi_par_i'*dF_damping;
matlabFunction(D_s, 'Vars', {q, dq, params, s}, "File", destination_folder_b + "/D_s");
D_q = str2func(sprintf("@(q, dq, params) myintegral(@(s) %s.D_s(q, dq, params, s)  , 0, params(1), 'ArrayValued', 1, 'AbsTol', %f, 'RelTol', %f)", className, IntAbsTol, IntRelTol));