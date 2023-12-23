function A = devectorize_3x3_matrix(a)
%DEVECTORIZE_3X3_MATRIX Transforms back the 6D vector a into a 3x3
%symmetrix matrix A. See also vectorize_3x3_matrix.m.
A = [a(1), a(4), a(5);
     a(4), a(2), a(6);
     a(5), a(6), a(3)];
end