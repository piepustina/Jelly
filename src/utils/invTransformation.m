function T_inv = invTransformation(T)
%INVTRANSFORMATION Computes the inverse of a transformation matrix
R       = T(1:3, 1:3);
t       = T(1:3, 4);
T_inv   = [R', -R'*t; zeros(1, 3), 1];
end

