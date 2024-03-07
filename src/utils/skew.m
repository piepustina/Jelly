function S = skew(X)
%SKEW Function that implements the skew symmetric operator
S = zeros(3, 3, "like", X);
S = [ 0   , -X(3),  X(2);
      X(3),  0   , -X(1);
     -X(2),  X(1),  0];
end

