function v = skew_inv(M)
%SKEW_INV Inverse skew symmetric operator for 3x3 matrices

arguments (Input)
    M (3, 3, :)
end

arguments (Output)
    v (3, :)
end

nM  = size(M, 3);
v   = [reshape(M(3, 2, 1:nM), 1, nM); ...
       reshape(M(1, 3, 1:nM), 1, nM) ; ...
       reshape(M(2, 1, 1:nM), 1, nM)];
end

