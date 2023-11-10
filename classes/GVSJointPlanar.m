classdef GVSJointPlanar < GVSJoint
    %GVSBODY Class representing a slender body modeled under the geometric
    %variable strain approach.
    
    %Specify number of DoFs
    properties(Constant)
        n    = 1;
    end
    
    %Implement the strain basis
    methods
        function obj = GVSJointPlanar(Parameters)
            obj      = obj@GVSJoint(GVSJointPlanar.n, Parameters);
        end
        function Phi = StrainBasis(obj, s)
            Phi = [0;-1/obj.RestLength;0;0;0;0];
        end
    end
end

