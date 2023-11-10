function adxi = adT(xi)
%AD Coadjoint operator
adxi = [skew(xi(1:3))   , skew(xi(4:6)); 
        zeros(3,3)      , skew(xi(1:3))];
end

