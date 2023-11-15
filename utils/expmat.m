function expmA = expmat(A)
%EXPMAT Matrix exponential for a skew symmetric matrix A
Theta = norm(skew_inv(A(1:3, 1:3)));

%Approximate theta if too small
if isnumeric(Theta)
    if Theta < 1e-5
        Theta = 1e-5;
    end
end

%Compute the exponential
expmA = eye(4) + A + (1-cos(Theta))/(Theta^2)*mpower(A, 2) + (Theta - sin(Theta))/(Theta^3)*mpower(A, 3);
end

