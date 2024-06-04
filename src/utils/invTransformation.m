function T_inv = invTransformation(T)
%INVTRANSFORMATION Computes the inverse of a transformation matrix or a 3D tensor of transformation matrices.

arguments (Input)
    T (4, 4, :)
end

arguments (Output)
    T_inv (4, 4, :)
end

% Number of transformation matrices
nT = size(T, 3);

% Extract the rotational part
R       = T(1:3, 1:3, 1:nT);
RT      = pagetranspose(R);
% Extract the translational part
t       = T(1:3, 4, 1:nT);
% Compute the inverse
T_inv   = [RT            , -pagemtimes(RT, t); ...
          zeros(1, 3, nT), ones(1, 1, nT)];
end

