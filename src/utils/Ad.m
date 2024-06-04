function Adg = Ad(g)
%AD Adjoint map of a transformation matrix g or a 3D tensor of transformation matrices g(:, :, i)

arguments (Input)
    g (4, 4, :)
end


arguments (Output)
    Adg (6, 6, :)
end

% Number of transformation matrices
nG      = size(g, 3);

% Extract the needed terms
Zeros   = zeros(3, 3, nG);
R       = g(1:3, 1:3, 1:nG);
r       = squeeze(g(1:3, 4, 1:nG));

% Compute the matrix adjoint representation
Adg     = [R                        , Zeros; 
           pagemtimes(skew(r), R)   , R];
end

