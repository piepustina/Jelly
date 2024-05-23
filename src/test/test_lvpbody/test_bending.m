% Test that the bending primitive is locally volume preserving
clear; clc; 

syms x [3, 1] real
syms q [1, 1] real
syms xi_A1(q) xi_A2(q) xi_A3(q)

syms xi_A2_v real

%%
alpha = ((1+xi_A1*x(2))-sqrt( (1+xi_A1*x(2))^2 -2*xi_A2*x(1)))/xi_A2;

dalpha_dx1 = simplify(diff(alpha, x(1)));
dalpha_dx2 = simplify(diff(alpha, x(2)));
hat_xi_A   = skew([xi_A1; xi_A2; xi_A3]);

I3         = eye(3);

J_f        = [I3(:, 1)*dalpha_dx1, I3(:, 1)*dalpha_dx2 + I3(:, 2), I3(:, 3) + hat_xi_A*[alpha; x(2); 0]];

% Compute the determinant of J_f, which must equal 1
simplify(det(J_f))


%%
d2 = xi_A1/xi_A2-xi_A1/xi_A2*(1+xi_A1*x(2))/sqrt((1+xi_A1*x(2))^2 - 2*xi_A2*x(1));





