function skewA = pageskew(A, transpose)
%#codegen

%PAGESKEW Applies the skew symmetrix operator to the columns of the matrix
%A \in R^{3 x n} and returns a 3D matrix of dimension 3 x 3 x n.

if nargin == 1
    transpose = 0;
end

n = size(A, 2);
skewA = zeros([3, 3, n], class(A));

for i = 1:n
    skewA(:, :, i) = skew(A(:, i));
    if transpose
        skewA(:, :, i) = skewA(:, :, i)';
    end
end

