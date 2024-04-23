classdef PGC3D < GVSBody
    %Class modeling a 3D piecewise Gaussian curvature (PGC) body with elongation.
    
    properties
        n    = 5;
    end
    
    methods
        function obj = PGC3D(Parameters)
            %Construct a 3D PGC body with elongation.
            %
            %Args:
            %   Parameters ([double], [sym]): Parameters of the body, specified as for :class:`GVSBody`
            obj      = obj@GVSBody(Parameters);
        end
        
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

