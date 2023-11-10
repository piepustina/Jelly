classdef GVSJoint3D < GVSJoint
    %GVSBODY Class representing a slender body modeled under the geometric
    %variable strain approach.
    
    %Specify number of DoFs
    properties(Constant)
        n    = 3;
    end
    
    %Implement the strain basis
    methods
        function obj = GVSJoint3D(Parameters)
            obj      = obj@GVSJoint(GVSJoint3D.n, Parameters);
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

