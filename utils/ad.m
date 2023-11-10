function adxi = ad(xi)
%AD Adjoint operator
adxi = [skew(xi(1:3)), zeros(3,3); skew(xi(4:6)), skew(xi(1:3))];
end

