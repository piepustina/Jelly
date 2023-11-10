function [C] = sympagemtimes(A, B)
%SYMPAGEMTIMES Implementation of the mpagemtimes function of Maltab for
%matrices containing symbolic elements.

%No check that the dimensions are correct.

[mA, nA, oA] = size(A);
[mB, nB, oB] = size(B);

%B is a 1/2D matrix
if oB == 1
    C = zeros([mA, nB, oA], class(A));
    for i = 1:oA
        C(:, :, i) = A(:, :, i)*B;
    end
    return;
end
%A is a 1/2D matrix
if oA == 1
    C = zeros([mA, nB, oB], class(B));
    for i = 1:oB
        C(:, :, i) = A*B(:, :, i);
    end
    return;
end

%We must have A and B with the same dimension
C = zeros([mA, nA, oA], class(A));
for i = 1:oA
    C(:, :, i) = A(:, :, i)*B(:, :, i);
end

end

