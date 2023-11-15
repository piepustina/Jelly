classdef PAC2D < GVSBody
    %Class representing a planar slender body modeled under the piecewise affine
    %curvature (PAC) hypothesis without elongation. 
    
    %Number of DoFs
    properties(Constant)
        n    = 2;
    end
    
    methods
        %Class constructor
        function obj = PAC2D(Parameters)
            obj      = obj@GVSBody(PAC2D.n, Parameters);
        end

        %Strain basis
        function Phi = StrainBasis(obj, s)
            Phi = [                0,                 0;
                   -1/obj.RestLength, -s/obj.RestLength;
                                   0,                 0;
                                   0,                 0;
                                   0,                 0;
                                   0,                 0];
        end
    end

end

