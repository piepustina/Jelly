function F = int_approx(f, x, intv, x_gauss, w_gauss)
%INT_APPROX Approximate a given integral using the Gaussian quadrature rule
% 

%Get the number of Gaussian quadrature points
n_Gauss = length(x_gauss);

%Get the extremal of integration
a = intv(1);
b = intv(2);

%Compute the approximated integral
F = 0;
for i = 1:n_Gauss
    F = F + w_gauss(i)*subs(f, x, (b - a)/2*x_gauss(i) + (a + b)/2);
end

F = (b - a)/2*F;

end

