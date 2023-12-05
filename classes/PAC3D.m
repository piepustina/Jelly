classdef PAC3D < GVSBody
    %Class representing a slender body modeled under the piecewise affine
    %curvature (PAC) hypothesis with elongation. 
    
    %Number of DoFs
    properties(Constant)
        n    = 5;
    end
    
    methods
        %Class constructor
        function obj = PAC3D(Parameters)
            obj      = obj@GVSBody(PAC3D.n, Parameters);
        end

        %Strain basis
        function Phi = StrainBasis(obj, s)
            L0 = obj.RestLength;
            Phi = [-1/L0, -s/L0,    0,    0, 0;
                       0,     0, 1/L0, s/L0, 0;
                       0,     0,    0,    0, 0;
                       0,     0,    0,    0, 0;
                       0,     0,    0,    0, 0;
                       0,     0,    0,    0, 1/L0];
        end
    end

end

