%Compute the model of rotational body (weird... to be changed).
%Make sure that the toolbox has been initialized

clear, clc;

%Define the configuration and its time derivatives
syms theta dtheta ddtheta real

%Define the parameters of the body:
% DH = DH parameters as  alpha|a|d|theta
% m = mass
% pcom = position of the CoM in the body frame
% I = inertia in the body frame

syms alpha_ a_ d_ m_ s real %The variable s is not used but it is required by the class.
pcom_ = sym('pcom_', [3, 1]);
I_ = sym('I_', [3,3]);
I_(2,1) = I_(1, 2);
I_(3,1) = I_(1, 3);
I_(3,2) = I_(2, 3);

%Parameters of the body
params = [alpha_; a_; d_; m_; pcom_; vectorize_3x3_matrix(I_)];

%Mass of the body
m = m_;

%Configuration variables
q = theta;
dq = dtheta;
ddq = ddtheta;


R_s = [cos(theta), -cos(alpha_)*sin(theta),  sin(alpha_)*sin(theta);
       sin(theta),  cos(alpha_)*cos(theta), -sin(alpha_)*cos(theta);
                0,             sin(alpha_),              cos(alpha_)];
t_s = [a_*cos(theta);
       a_*sin(theta);
       d_];

T_s  = [R_s, t_s
     zeros(1, 3), 1];

T  = T_s;

%Compute the strain field
xi_hat = zeros(4, 4, 'sym');
%Extract the strain as a 6D vector in angular and linear part
xi = [skew_inv(xi_hat(1:3, 1:3)); xi_hat(1:3, 4)];

T_inv = simplify(inv(T));

%Take the rotation matrix and the position vector
R_i_1_i = T(1:3, 1:3);
t_i_1_i = T(1:3, 4);

prev_p_i = T_s(1:3, 4);

%Rotate the position vector in {S_i}
p_i    = simplify(T_inv*[prev_p_i;1]);
p_i    = p_i(1:3);
p_comi = pcom_;

%Compute the relative position vector and its derivatives (which should be zero)
r_i    = simplify(p_i - p_comi);
dr_i   = simplify(matrix_derivative(r_i, q, dq));
ddr_i  = simplify(matrix_derivative(dr_i, [q;dq], [dq;ddq]));

%% Compute all the required terms
destination_folder_b = "./classes/@Rotational_Body";
destination_folder_j = "./classes/@Rotational_Joint";


%% Compute the kinematic data
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
v_par_i   = simplify(R_i_1_i'*jacobian(v_i_1_i, dq));
%Partial angular velocity
omega_par_i   = simplify(R_i_1_i'*jacobian(omega_i_1_i, dq));
%Strain field
xi_i = xi;
%Strain field in the reference configuration
xi_i_star = subs(xi_i, q, zeros(size(q)));
%Partial strain
xi_par_i  = jacobian(xi, q);

%%
I = I_;
dI = zeros(3, 3, 'sym');
J = zeros(3, 3, 'sym');
int_dr_i = zeros(3, 1, 'sym');
int_ddr_i = zeros(3, 1, 'sym');
int_r_i_X_dr_i = zeros(3, 1, 'sym');
int_r_i_X_ddr_i = zeros(3, 1, 'sym');
int_dr_i_X_pv_r = zeros(size(jacobian(dr_i, dq)'*skew(dr_i)'), 'sym');
int_pv_r_O_dd_r = zeros(size(jacobian(dr_i, dq)'*ddr_i), 'sym');
int_dr_i_O_dr_i = zeros([1,1], 'sym');
grad_int_dr_i = zeros(size(jacobian(dr_i, dq)'), 'sym');
grad_int_r_i_X_dr_i = zeros(size(jacobian(cross(r_i, dr_i), dq)'), 'sym');
grad_J = zeros(3, 3, length(q), 'sym');
grad_v_com_i = jacobian(v_i_comi, dq)';
K_q     = zeros(1, 1, 'sym');
D_q     = zeros(1, 1, 'sym');

%% Export all the functions
disp("Exporting functions...");
cF = pwd;
cd(destination_folder_j);

matlabFunction(T                   , 'Vars', {q, params}, 'File', "T");
matlabFunction(T_s                 , 'Vars', {q, params, s}, 'File', "T_s");
matlabFunction(v_i_1_i             , 'Vars', {q, dq, params}, 'File', "v_rel");
matlabFunction(omega_i_1_i         , 'Vars', {q, dq, params}, 'File', "omega_rel");
matlabFunction(dv_i_1_i            , 'Vars', {q, dq, ddq, params}, 'File', "a_rel");
matlabFunction(domega_i_1_i        , 'Vars', {q, dq, ddq, params}, 'File', "domega_rel");
matlabFunction(v_par_i             , 'Vars', {q, params}, 'File', "v_par");
matlabFunction(omega_par_i         , 'Vars', {q, params}, 'File', "omega_par");

cd(cF)
cd(destination_folder_b);

matlabFunction(p_comi              , 'Vars', {q, params}, 'File',"p_com");
matlabFunction(v_i_comi            , 'Vars', {q, dq, params}, 'File', "v_com_rel");
matlabFunction(dv_i_comi           , 'Vars', {q, dq, ddq, params}, 'File', "a_com_rel");
matlabFunction(I                   , 'Vars', {q, params}, 'File', "I");
matlabFunction(m                   , 'Vars', {params}, 'File', "m");
matlabFunction(dI                  , 'Vars', {q, dq, params}, 'File', "dI");
matlabFunction(J                   , 'Vars', {q, dq, params}, 'File', "J");
matlabFunction(int_dr_i            , 'Vars', {q, dq, params}, 'File', "int_dr");
matlabFunction(int_ddr_i           , 'Vars', {q, dq, ddq,  params}, 'File', "int_ddr");
matlabFunction(int_r_i_X_dr_i      , 'Vars', {q, dq,  params}, 'File', "int_r_X_dr");
matlabFunction(int_r_i_X_ddr_i     , 'Vars', {q, dq, ddq,  params}, 'File', "int_r_X_ddr");
matlabFunction(int_dr_i_X_pv_r     , 'Vars', {q, dq,  params}, 'File', "int_dr_X_pv_r");
matlabFunction(int_pv_r_O_dd_r     , 'Vars', {q, dq, ddq,  params}, 'File', "int_pv_r_O_dd_r");
matlabFunction(int_dr_i_O_dr_i     , 'Vars', {q, dq,  params}, 'File', "int_dr_O_dr");
matlabFunction(grad_int_dr_i       , 'Vars', {q,  params}, 'File', "grad_int_dr");
matlabFunction(grad_int_r_i_X_dr_i , 'Vars', {q,  params}, 'File', "grad_int_r_X_dr");
matlabFunction(grad_J              , 'Vars', {q,  params}, 'File', "grad_J");
matlabFunction(grad_v_com_i        , 'Vars', {q,  params}, 'File', "grad_v_com");
matlabFunction(xi_i                , 'Vars', {q, params, s}, 'File', 'xi');
matlabFunction(K_q                 , 'Vars', {q, params}, 'File', 'K');
matlabFunction(D_q                 , 'Vars', {q, dq, params}, 'File', 'D');
cd(cF)