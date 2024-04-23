function adxi = ad(xi)
%AD Adjoint operator
%adxi = zeros(6, 6, "like", xi);
% adxi = [skew(xi(1:3)), zeros(3,3); ...
%         skew(xi(4:6)), skew(xi(1:3))];

arguments (Input)
    xi (6, :) double
end

arguments (Output)
    adxi (6, 6, :) double
end

% Vectorized version
nxi  = size(xi, 2);
adxi = zeros(6, 6, nxi, "like", xi);

% Assign the output
Sxi13 = skew(xi(1:3, 1:nxi));
adxi(1:3, 1:3, 1:nxi) = Sxi13;
adxi(4:6, 4:6, 1:nxi) = Sxi13;
adxi(4:6, 1:3, 1:nxi) = skew(xi(4:6, 1:nxi));
end

