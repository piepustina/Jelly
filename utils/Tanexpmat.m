function TA = Tanexpmat(A)
%TANMAT Compute the tangent operator to the exponential map of the skew
%symmetric matrix A

Theta = norm(A(1:3));

%Approximate theta if too small
if isnumeric(Theta)
    if Theta < 1e-6
        Theta = 1e-6;
    end
end

adA   = ad(A);
TA    = eye(6) + (4-4*cos(Theta)-Theta*sin(Theta))/(2*Theta^2)*adA+...
                 (4*Theta-5*sin(Theta)+Theta*cos(Theta))/(2*Theta^3)*mpower(adA, 2)+...
                 (2-2*cos(Theta)-Theta*sin(Theta))/(2*Theta^4)*mpower(adA, 3)+...
                 (2*Theta-3*sin(Theta)+Theta*cos(Theta))/(2*Theta^5)*mpower(adA, 4);

% %Numerical approximation using Gaussian quadrature rule
% N = 50;
% [xGauss, wGauss] = lgwt(N, 0, 1);
% TA = zeros(6, 6);
% for i = 1:N
%     TA = TA + wGauss(i)*expm(xGauss(i)*adA);
% end
end

