function Adg = invAd(g)
%AD Inverse Adjoint map
Zeros   = zeros(3, 3);
RT      = g(1:3, 1:3)';
r       = g(1:3, 4);
Adg     = [ RT          , Zeros; 
            -RT*skew(r) , RT];
end

