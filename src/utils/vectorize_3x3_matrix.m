function a = vectorize_3x3_matrix(A)
%VECTORIZE_3X3_MATRIX Vectorizes a 3x3 symmetric matrix into a 6d vector
a = [A(1, 1); A(2, 2); A(3, 3); A(1, 2); A(1, 3); A(2, 3)];
end

