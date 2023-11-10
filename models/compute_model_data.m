%%
%--------------------------------------------------------------------------
% Compute relative velocities and accelerations of the planar PCC joint----
%--------------------------------------------------------------------------
disp("Computing kinematic data...");
%Relative linear velocity
v_i_1_i = matrix_derivative(t_i_1_i, q, dq);
v_i_1_i = simplify(v_i_1_i);
%Relative angular velocity
S_omega_i_1_i = matrix_derivative(R_i_1_i, q, dq)*R_i_1_i';
omega_i_1_i   = [S_omega_i_1_i(3,2);S_omega_i_1_i(1,3);S_omega_i_1_i(2,1)];
omega_i_1_i = simplify(omega_i_1_i);
%Relative linear acceleration
dv_i_1_i = matrix_derivative(v_i_1_i, [q;dq], [dq;ddq]);
%Relative angular acceleration
domega_i_1_i = matrix_derivative(omega_i_1_i, [q;dq], [dq;ddq]);
%Relative CoM velocity
v_i_comi = matrix_derivative(p_comi, q, dq);
%Relative CoM acceleration
dv_i_comi = matrix_derivative(v_i_comi, [q; dq], [dq; ddq]);
%Partial linear velocity
v_par_i   = R_i_1_i'*jacobian(v_i_1_i, dq);
%Partial angular velocity
omega_par_i   = R_i_1_i'*jacobian(omega_i_1_i, dq);
%Strain field
xi_i = xi;
%Strain field in the reference configuration
xi_i_star = xi_i;
for i = 1:length(q)
    xi_i_star = limit(xi_i_star, q(i), 0);
end
%Partial strain
xi_par_i  = jacobian(xi, q);
%%
%--------------------------------------------------------------------------
% Compute all the needed integrals ----------------------------------------
%--------------------------------------------------------------------------
disp("Computing integrals...");
disp("Computing I");
I = int(simplify(skew(r_i)'*skew(r_i)*rho_L),...
        s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
I = simplify(I);
%--------------------------------------------------------------------------
disp("Computing dI");
% dI = simplify(matrix_derivative(I, q, dq));
% dI = int(matrix_derivative(skew(r_i)'*skew(r_i)*rho_L, q, dq), ...
%         s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
% dI = simplify(dI);
dI = zeros(3, 3, 'sym');


%--------------------------------------------------------------------------
disp("Computing J");
J = int(simplify(skew(dr_i)'*skew(r_i)*rho_L + skew(r_i)'*skew(dr_i)*rho_L),...
        s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
J = simplify(J);

%--------------------------------------------------------------------------
disp("Computing int_dr_i");
% int_dr_i = int(dr_i*rho_L ,...
%                     s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
% int_dr_i = simplify(int_dr_i);

int_dr_i = zeros(3, 1, 'sym');

%--------------------------------------------------------------------------
disp("Computing int_ddr_i");

% int_ddr_i = int(ddr_i*rho_L,...
%                 s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
% int_ddr_i = simplify(int_ddr_i);

int_ddr_i = zeros(3, 1, 'sym');

%--------------------------------------------------------------------------
disp("Computing int_r_i_X_dr_i");

int_r_i_X_dr_i = int(simplify(cross(r_i, dr_i)*rho_L),...
                        s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
int_r_i_X_dr_i = simplify(int_r_i_X_dr_i);
%--------------------------------------------------------------------------
disp("Computing int_r_i_X_ddr_i");
int_r_i_X_ddr_i = int(simplify(cross(r_i, ddr_i)*rho_L),...
                        s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
int_r_i_X_ddr_i = simplify(int_r_i_X_ddr_i);
%--------------------------------------------------------------------------
disp("Computing int_dr_i_X_pv_r");
int_dr_i_X_pv_r = int(jacobian(dr_i, dq)'*skew(dr_i)'*rho_L, ...
                        s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
int_dr_i_X_pv_r = simplify(int_dr_i_X_pv_r);
%--------------------------------------------------------------------------
disp("Computing int_pv_r_O_dd_r");

int_pv_r_O_dd_r = int(jacobian(dr_i, dq)'*ddr_i*rho_L, ...
                        s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
int_pv_r_O_dd_r = simplify(int_pv_r_O_dd_r);
%--------------------------------------------------------------------------
disp("Computing int_dr_i_O_dr_i");

% int_dr_i_O_dr_i = int(dr_i'*dr_i*rho_L, ...
%                         s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
% int_dr_i_O_dr_i = simplify(int_dr_i_O_dr_i);

int_dr_i_O_dr_i = zeros([1,1], 'sym');

%%
%--------------------------------------------------------------------------
% Compute all the needed gradients ----------------------------------------
%--------------------------------------------------------------------------
disp("Computing gradients...");

% Since all these depend on dq, we can bring the differential operator
% inside the integrals!
% grad_int_dr_i = jacobian(int_dr_i, dq)';
% grad_int_r_i_X_dr_i = jacobian(int_r_i_X_dr_i, dq)';
% 
% grad_J = zeros(3, 3, length(q), 'sym');
% for i = 1:length(q)
%     grad_J(:, :, i) = diff(J, dq(i));
% end
% 
% grad_int_dr_i = int(jacobian(dr_i*rho_L, dq)' ,...
%                     s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);

grad_int_dr_i = zeros(size(jacobian(dr_i*rho_L, dq)'), 'sym');

grad_int_r_i_X_dr_i = int(simplify(jacobian(cross(r_i, dr_i)*rho_L, dq)'),...
                        s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
grad_J = zeros(3, 3, length(q), 'sym');
for i = 1:length(q)
    grad_J(:, :, i) = int(simplify(diff(skew(dr_i)'*skew(r_i)*rho_L + skew(r_i)'*skew(dr_i)*rho_L, dq(i))),...
                            s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
end

grad_v_com_i = jacobian(v_i_comi, dq)';

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

%Body damping matrices
D_l = K_l*Eta;
D_a = K_a*Eta;
%Body damping matrix
D_b = blkdiag(D_a, D_l);

%Body damping force
dF_damping = D_b*matrix_derivative(xi_i, q, dq);

%Computation of the generalized forces by projecting the forces into the
%strain space
K_s = xi_par_i'*dF_elastica;
K_q = int(K_s, s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);
D_s   = xi_par_i'*dF_damping;
D_q = int(D_s, s, [0, L], 'IgnoreSpecialCases', true, 'Hold', hold_integrals);