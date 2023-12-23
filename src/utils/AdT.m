function Adg = AdT(g)
%AD Coadjoint map
Zeros   = zeros(3, 3);
R       = g(1:3, 1:3);
r       = g(1:3, 4);
Adg     = [ R           , skew(r)*R; 
            Zeros       , R];
end

