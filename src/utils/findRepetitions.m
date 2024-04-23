function [Y, C] = findRepetitions(X)
%FINDREPETITIONS Given a vector X finds the number of consecutive repetitions. The function also returns
% a compact representation C of X such that X = repelem(C, Y)

arguments (Input)
    X (1, :) double
end

arguments (Output)
    Y (1, :) double
    C (1, :) double
end

d = [true, diff(X) ~= 0, true];  % TRUE if values change
Y = diff(find(d));               % Number of repetitions for each element

lX = length(X);
C = X(d(1:lX));
end

