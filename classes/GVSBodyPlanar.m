classdef GVSBodyPlanar < GVSBody
    %GVSBODY Class representing a slender planar body modeled under PCC
    %hypothesis.
    
    %Number of DoFs
    properties(Constant)
        n    = 1;
    end
    
    %Implement the strain basis
    methods
        function obj = GVSBodyPlanar(Parameters)
            obj      = obj@GVSBody(GVSBodyPlanar.n, Parameters);
        end

        function Phi = StrainBasis(obj, s)
            Phi = [0;-1/obj.RestLength;0;0;0;0];
        end
    end

end

