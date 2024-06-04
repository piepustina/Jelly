function v = skew4_inv(M)
%SKEW_INV Inverse skew symmetric operator for 4x4 matrices

arguments (Input)
    M (4, 4, :)
end

arguments (Output)
    v (6, :)
end

% Compute the number of matrices
nM = size(M, 3);

% Asssign the output
v = [M(3, 2, 1:nM); ...
     M(1, 3, 1:nM); ...
     M(2, 1, 1:nM); ...
     M(1:3, 4, 1:nM)];

end

