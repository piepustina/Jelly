%Run this script to compute the model of a PCC planar slender soft body and
%of its joint. 
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

T_s = [cos(s/L*theta), -sin(s/L*theta), 0, L*sin(s/L*theta)/theta    ;
     sin(s/L*theta),  cos(s/L*theta), 0,   L*(1-cos(s/L*theta))/theta;
              0,           0, 1,                      0;
              0,           0, 0,                      1];

T_s = [cos(s/L*theta), 0, -sin(s/L*theta), -L*(1-cos(s/L*theta))/theta;
       0             , 1,               0, 0;
       sin(s/L*theta), 0,  cos(s/L*theta), L*sin(s/L*theta)/theta    ;
              0,       0, 0,                      1];

T = subs(T_s, s, L);

            
T_inv = simplify(invTransformation(T));            
disp("The transformation matrix is:");
pretty(T)

%Take the rotation matrix and the position vector
R_i_1_i = T(1:3, 1:3);
t_i_1_i = T(1:3, 4);


%Compute the strain field
xi_hat = simplify(T_s\diff(T_s, s));

%Extract the strain as a 6D vector in angular and linear part
xi = [skew_inv(xi_hat(1:3, 1:3)); xi_hat(1:3, 4)];


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

destination_folder_j = "./oldClasses/@PCCPlanar_Joint_no_elongation";
destination_folder_b = "./oldClasses/@PCCPlanar_Body_no_elongation";

compute_model_data;

toc

disp("Export the functions...");
export_functions;