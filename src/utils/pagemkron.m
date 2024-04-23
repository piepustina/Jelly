%#codegen
function K = pagemkron(A, B)
%PAGEMKRON Compute the kronecker product between the pages of the 3D tensors A and B. 

% Get the sizes of the inputs
ma = size(A, 1);
na = size(A, 2);
oa = size(A, 3);
mb = size(B, 1);
nb = size(B, 2);
ob = size(B, 3);

% Check that the third dimension of A and B is the same
if oa ~= ob
    error("The third size of A and B must be equal.");
end

% Compute the Kronecker product between the pages
A = reshape(A, 1, ma, 1, na, oa);
B = reshape(B, mb, 1, nb, 1, ob);
K = reshape(A.*B, ma*mb, na*nb, oa);
end

