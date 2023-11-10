%Run this script to compute the model of a PCC planar slender soft body and
%of its joint. This script uses the gamma function representation of the
%sinc function. In particular, we use the identities
% sinc(x) = sin(x)/x = 1/(gamma(1-x/pi)*gamma(1+x/pi))
% (1-cos(x))/x = sin(x/2)*sin(x/2)/(x/2) = sin(x/2)*sinc(x/2)
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

%%
%Set to true to leave the integrals in implicit form
hold_integrals = false;
%Override the simplify function to speedup computation
%simplify = @(varargin) simplify(varargin{:});
simplify = @(varargin) varargin{1};
%Override the integral operator
% int      = @(varargin) int(varargin{:});
%Number of Gaussian points to approximate the integrals
n_Gauss = 3;
syms x real;
%Find the roots of the legendre polynomial
legendre_poly  = legendreP(n_Gauss, x);
dlegendre_poly = diff(legendre_poly, x);
x_gauss = vpasolve(legendre_poly == 0);
w_gauss = vpa(2./((1-x_gauss.^2).*subs(dlegendre_poly, x, x_gauss).^2), 4);
int      = @(varargin) int_approx(varargin{1}, varargin{2}, varargin{3}, x_gauss, w_gauss);



%%            
T_inv = simplify(invTransformation(T));            
disp("The transformation matrix is:");
pretty(T)

%Take the rotation matrix and the position vector
R_i_1_i = T(1:3, 1:3);
t_i_1_i = T(1:3, 4);


%Compute the strain field
xi_hat = simplify(T_s_strain\diff(T_s_strain, s));

%Extract the strain as a 6D vector in angular and linear part
xi = [skew_inv(xi_hat(1:3, 1:3)); xi_hat(1:3, 4)];

prev_p_i = T_s(1:3, 4);

%Compute the CoM in {S_{i-1}}
prev_p_comi = (1/m)*int(prev_p_i*rho_L, s, [0, L]);


%prev_p_comi = simplify(prev_p_comi);
%prev_r_i    = simplify(prev_p_i - prev_p_comi);

%Rotate the vectors in {S_i}
p_comi = T_inv*[prev_p_comi;1];
p_i    = T_inv*[prev_p_i;1];

%p_comi = simplify(p_comi(1:3));
%p_i    = simplify(p_i(1:3));
p_comi = p_comi(1:3);
p_i    = p_i(1:3);

%Compute the relative position vector and its derivatives
r_i    = p_i - p_comi;
dr_i   = matrix_derivative(r_i, q, dq);
ddr_i  = matrix_derivative(dr_i, [q;dq], [dq;ddq]);

%assert(~all(simplify(int(r_i*rho_L, s, [0, L]))))

%%

%Compute all the required terms
destination_folder_j = "./classes/@PCCPlanar_Joint_no_elongation_Gamma";
destination_folder_b = "./classes/@PCCPlanar_Body_no_elongation_Gamma";

compute_model_data;

toc

disp("Export the functions...");
export_functions;