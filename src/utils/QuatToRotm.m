%#codegen
function R = QuatToRotm(q)
%QUATTOROTM Convert a quaternion to a rotation matrix

arguments (Input)
    q (4, :) double
end

arguments (Output)
    R (3, 3, :) double
end

% Preallocate the output
nR = size(q, 2);
R  = zeros(3, 3, nR);

qEta    = q(1, 1:nR);
qX      = q(2, 1:nR);
qY      = q(3, 1:nR);
qZ      = q(4, 1:nR);

% Build the diagonal
R(1, 1, 1:nR) = 2.*(qEta.^2 + qX.^2) - 1;
R(2, 2, 1:nR) = 2.*(qEta.^2 + qY.^2) - 1;
R(3, 3, 1:nR) = 2.*(qEta.^2 + qZ.^2) - 1;

% Build the upper diagonal part
R(1, 2, 1:nR) = 2.*(qX.*qY - qEta.*qZ);
R(1, 3, 1:nR) = 2.*(qX.*qZ + qEta.*qY);
R(2, 3, 1:nR) = 2.*(qY.*qZ - qEta.*qX);

% Build the lower diagonal part
R(2, 1, 1:nR) = 2.*(qX.*qY + qEta.*qZ);
R(3, 1, 1:nR) = 2.*(qX.*qZ - qEta.*qY);
R(3, 2, 1:nR) = 2.*(qY.*qZ + qEta.*qX);

end

