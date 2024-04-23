function expmA = expmat(A)
%EXPMAT Matrix exponential for a skew symmetric matrix A

% % Nonvectiruzed version
% Theta = norm(skew_inv(A(1:3, 1:3)));
% 
% %Approximate theta if too small
% if isnumeric(Theta)
%     if Theta < 1e-6
%         Theta = 1e-6;
%     end
% end
% 
% %Compute the exponential
% expmA = eye(4) + A + (1-cos(Theta))/(Theta^2)*mpower(A, 2) + (Theta - sin(Theta))/(Theta^3)*mpower(A, 3);

arguments (Input)
    A (4, 4, :)
end

arguments (Output)
    expmA (4, 4, :)
end

% Vectorized version
nA    = size(A, 3);
Theta = reshape(vecnorm(skew_inv(A(1:3, 1:3, 1:nA))), 1, 1, nA);

%Approximate theta if too small
if isnumeric(Theta)
    % Approximate Theta to compute the limit
    Theta(Theta < 1e-6) = 1e-6;
end

%Compute the exponential
A2 = pagemtimes(A, A);
A3 = pagemtimes(A2, A);

expmA = repmat(eye(4), 1, 1, nA) ...
      + A ...
      + pagemtimes((1-cos(Theta))./(Theta.^2), A2) ...
      + pagemtimes((Theta - sin(Theta))./(Theta.^3), A3);
end

