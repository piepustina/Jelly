%#codegen
function [q] = RotmToQuat(R)
%ROTTOQUAT Convert a rotation matrix to quaternion format

arguments (Input)
    R (3, 3, :) double
end

arguments (Output)
    q (4, :)
end

% Output preallocation
nR = size(R, 3);
q = zeros(4, nR);

% Compute the quaternions
R11        = R(1, 1, 1:nR);
R22        = R(2, 2, 1:nR);
R33        = R(3, 3, 1:nR);
q(1, 1:nR) = 0.5.*real(sqrt(1 + R11 + R22 + R33));
q(2, 1:nR) = 0.5.*sgn(R(3, 2, 1:nR) - R(2, 3, 1:nR)).*real(sqrt(R11 - R22 - R33 + 1));
q(3, 1:nR) = 0.5.*sgn(R(1, 3, 1:nR) - R(3, 1, 1:nR)).*real(sqrt(R22 - R11 - R33 + 1));
q(4, 1:nR) = 0.5.*sgn(R(2, 1, 1:nR) - R(1, 2, 1:nR)).*real(sqrt(R33 - R22 - R11 + 1));

end


