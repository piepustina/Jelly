function v = skew4_inv(M)
%SKEW_INV Inverse skew symmetric operator for 4x4 matrices
v = [M(3,2); M(1,3); M(2,1); M(1:3, 4)];
end

