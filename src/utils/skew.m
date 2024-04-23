%#codegen
function S = skew(X)
%SKEW Function that implements the skew symmetric operator

arguments (Input)
    X (3, :)
end

arguments (Output)
    S (3, 3, :)
end

% Output preallocation
nX = size(X, 2);
S = zeros(3, 3, nX, "like", X);

% Build the skew symmetric matrix in vectorized form
S(2, 1, 1:nX) =  X(3, 1:nX);
S(3, 1, 1:nX) = -X(2, 1:nX);
S(2, 3, 1:nX) = -X(1, 1:nX);
S(1, 2, 1:nX) = -S(2, 1, 1:nX);
S(1, 3, 1:nX) = -S(3, 1, 1:nX);
S(3, 2, 1:nX) = -S(2, 3, 1:nX);
end

