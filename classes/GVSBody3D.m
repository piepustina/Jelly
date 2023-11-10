classdef GVSBody3D < GVSBody
    %GVSBODY Class representing a slender planar body modeled under PCC
    %hypothesis.
    
    %Number of DoFs
    properties(Constant)
        n    = 3;
    end
    
    %Implement the strain basis
    methods
        function obj = GVSBody3D(Parameters)
            obj      = obj@GVSBody(GVSBody3D.n, Parameters);
        end

        function Phi = StrainBasis(obj, s)
            L0 = obj.RestLength;
            Phi = [   0,   -1/L0, 0;
                        1/L0,       0, 0;
                           0,       0, 0;
                           0,       0, 0;
                           0,       0, 0;
                           0,       0, 1/L0];
        end
    end

end

