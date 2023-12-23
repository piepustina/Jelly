%Compute the model of a PCC segment with cylindrical shape.
%All the inertial parameters, e.g., mass density, Young's modulus, etc...,
%are assumed constant.
%Make sure that the toolbox has been initialized

clear, clc;

%Define the configuration and its time derivatives
syms theta phi deltaL dtheta dphi ddeltaL ddtheta ddphi dddeltaL real

%Define the parameters of the body assumed linear isotropic:
% L_0 = Rest length
% rho = mass density of the body
% R   = cross sectional radius
% E = Young's modulus
% Poi = Poisson ration
% Eta = Material damping

syms L_0 R rho E Poi Eta real

%Reference length used to compute the integral
L = L_0;
%Mass density along the curvilinar abscissa
rho_L = rho*pi*R^2;

%Mass of the segment
m = rho_L*L_0;

q   = [theta; phi; deltaL];
dq  = [dtheta; dphi; ddeltaL];
ddq = [ddtheta; ddphi; dddeltaL];

params = [L_0; R; rho; E; Poi; Eta];

%Define the variables of integration
syms s real

%Delta parametrization
Delta_X = q(1);
Delta_Y = q(2);
Delta_L = q(3);
Delta   = sqrt(Delta_X^2 + Delta_Y^2);

R_s = [ 1 + (Delta_X/Delta)^2*(cos(s/L*Delta) - 1)    , (Delta_X*Delta_Y)/Delta^2*(cos(s/L*Delta) - 1), Delta_X/Delta*sin(s/L*Delta);
        (Delta_X*Delta_Y)/Delta^2*(cos(s/L*Delta) - 1), 1 + (Delta_Y/Delta)^2*(cos(s/L*Delta) - 1)    , Delta_Y/Delta*sin(s/L*Delta); 
        -Delta_X/Delta*sin(s/L*Delta)                  , -Delta_Y/Delta*sin(s/L*Delta)                  , cos(s/L*Delta)];
t_s = ((L_0 + Delta_L)/Delta^2)*[Delta_X*(1-cos(s/L*Delta)); Delta_Y*(1 - cos(s/L*Delta)); Delta*sin(s/L*Delta)];

%
T_s  = [R_s, t_s
     zeros(1, 3), 1];

T  = subs(T_s, s, L);

%Compute the strain field
xi_hat = simplify(T_s\diff(T_s, s));

%Extract the strain as a 6D vector in angular and linear part
xi = [skew_inv(xi_hat(1:3, 1:3)); xi_hat(1:3, 4)];

%%

T_inv = simplify(inv(T));            
disp("The transformation matrix is:");
pretty(T)

%Take the rotation matrix and the position vector
R_i_1_i = T(1:3, 1:3);
t_i_1_i = T(1:3, 4);

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

r_i = simplify(r_i);
dr_i = simplify(dr_i);
ddr_i = simplify(ddr_i);

assert(~all(simplify(int(r_i*rho_L, s, [0, L]))))

%Set to true to leave the integrals in implicit form
hold_integrals = false;
%Override the simplify function to speedup computation
%simplify = @(varargin) simplify(varargin{:});
simplify = @(varargin) varargin{1};
%Override the integral operator
% int      = @(varargin) int(varargin{:});
%Number of Gaussian points to approximate the integrals
n_Gauss = 5;
syms x real;
%Find the roots of the legendre polynomial
legendre_poly  = legendreP(n_Gauss, x);
dlegendre_poly = diff(legendre_poly, x);
x_gauss = vpasolve(legendre_poly == 0);
w_gauss = vpa(2./((1-x_gauss.^2).*subs(dlegendre_poly, x, x_gauss).^2), 4);
int      = @(varargin) int_approx(varargin{1}, varargin{2}, varargin{3}, x_gauss, w_gauss);

%Compute all the required terms
destination_folder_b = "./oldClasses/@PCCCylindricalDelta_Body";
destination_folder_j = "./oldClasses/@PCCCylindricalDelta_Joint";

compute_model_data;

disp("Export the functions...");
export_functions;