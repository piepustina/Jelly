function Adg = invAd(g)
%AD Inverse Adjoint map
% Non vectorized version
% Adg = zeros(6, 6, "like", g);
% Zeros   = zeros(3, 3);
% RT      = g(1:3, 1:3)';
% r       = g(1:3, 4);
% Adg     = [ RT          , Zeros; 
%             -RT*skew(r) , RT];

% Vectorized version

arguments (Input)
    g (4, 4, :)
end

arguments (Output)
    Adg (6, 6, :)
end

% Output preallocation
ng      = size(g, 3);
Adg     = zeros(6, 6, ng, "like", g);

% Store useful quantities
RT                  = pagetranspose(g(1:3, 1:3, 1:ng));
r                   = squeeze(g(1:3, 4, 1:ng));
% Assign the output
Adg(1:3, 1:3, 1:ng) = RT;
Adg(4:6, 4:6, 1:ng) = RT;
Adg(4:6, 1:3, 1:ng) = -pagemtimes(RT, skew(r));
end

