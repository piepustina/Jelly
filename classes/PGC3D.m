classdef PGC3D < GVSBody
    %GVSBODY Class representing a slender 3D body modeled under the Piecewise Gaussian Curvature (PGC)
    %hypothesis.
    
    %Number of DoFs
    properties(Constant)
        n    = 5;
    end
    
    methods
        %Class constructor
        function obj = PGC3D(Parameters)
            obj      = obj@GVSBody(PGC3D.n, Parameters);
        end
        
        %Strain basis
        function Phi = StrainBasis(obj, s)
            L0 = obj.RestLength;
            Phi = [-1/L0, -exp(-(s-L0/2)^2)/L0,    0,                 0, 0;
                       0,                  0, 1/L0, exp(-(s-L0/2)^2)/L0, 0;
                       0,                  0,    0,                 0, 0;
                       0,                  0,    0,                 0, 0;
                       0,                  0,    0,                 0, 0;
                       0,                  0,    0,                 0, 1/L0];
        end
    end

end

