function dM = matrix_derivative(M, var_, dvar_)
%DMDT Compute the time derivative of the symbolic matrix M(q) with dq =
%diff(q)/t
[m,n] = size(M);
dM = zeros([m,n], 'sym');

for i = 1:m
    for j = 1:n
        dM(i, j) = jacobian(M(i,j), var_)*dvar_;
    end
end

%dM = simplify(dM);
end

% %Second implementation (slower)
% function dM = matrix_derivative(M, var_, dvar_)
% %DMDT Compute the time derivative of the symbolic matrix M(q) with dq =
% %diff(q)/t
% [m,n] = size(M);
% dM = zeros([m,n], 'sym');
% 
% for i = 1:m
%     for j = 1:n
%         dM(i, j) = jacobian(M(i,j), var_)*dvar_;
%     end
% end
% 
% dM = simplify(dM);
% end
% 
