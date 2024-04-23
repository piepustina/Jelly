function TA = Tanexpmat(A)
%TANMAT Compute the tangent operator to the exponential map of the skew
%symmetric matrix A

% Non vectorized version
% Theta = norm(A(1:3));
% 
% %Approximate theta if too small
% if isnumeric(Theta)
%     if Theta < 1e-6
%         Theta = 1e-6;
%     end
% end
% 
% adA   = ad(A);
% TA    = eye(6) + (4-4*cos(Theta)-Theta*sin(Theta))/(2*Theta^2)*adA+...
%                  (4*Theta-5*sin(Theta)+Theta*cos(Theta))/(2*Theta^3)*mpower(adA, 2)+...
%                  (2-2*cos(Theta)-Theta*sin(Theta))/(2*Theta^4)*mpower(adA, 3)+...
%                  (2*Theta-3*sin(Theta)+Theta*cos(Theta))/(2*Theta^5)*mpower(adA, 4);

% %Numerical approximation using Gaussian quadrature rule
% N = 50;
% [xGauss, wGauss] = lgwt(N, 0, 1);
% TA = zeros(6, 6);
% for i = 1:N
%     TA = TA + wGauss(i)*expm(xGauss(i)*adA);
% end

% Vectorized version 
arguments (Input)
    A (6, :) double
end

arguments (Output)
    TA (6, 6, :) double
end

% Compute the third dimension of the output
nA = size(A, 2);

% Compute the angle
Theta = reshape(vecnorm(A(1:3, :)), 1, 1, []);

%Approximate theta if too small
if isnumeric(Theta)
    % Set Theta to a threshold for the numerical computations
    Theta(Theta < 1e-6) = 1e-6;
end

% Compute useful quantities
adA     = ad(A);
adA2    = pagemtimes(adA, adA);
adA3    = pagemtimes(adA2, adA);
adA4    = pagemtimes(adA3, adA);
cTheta  = cos(Theta);
sTheta  = sin(Theta);


TA    = repmat(eye(6), 1, 1, nA) ...
      + pagemtimes((4-4*cTheta-Theta.*sTheta)./(2*Theta.^2), adA) ...
      + pagemtimes((4*Theta-5*sTheta+Theta.*cTheta)./(2*Theta.^3), adA2) ...
      + pagemtimes((2-2*cTheta-Theta.*sTheta)./(2*Theta.^4), adA3) ...
      + pagemtimes((2*Theta-3*sTheta+Theta.*cTheta)./(2*Theta.^5), adA4);
end

