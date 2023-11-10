function Adg = Ad(g)
%AD Adjoint map
Zeros   = zeros(3, 3);
R       = g(1:3, 1:3);
r       = g(1:3, 4);
Adg     = [ R           , Zeros; 
            skew(r)*R   , R];
end

