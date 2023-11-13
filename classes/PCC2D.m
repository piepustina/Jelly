classdef PCC2D < GVSBody
    %GVSBODY Class representing a slender 2D body modeled under the PCC
    %hypothesis. The configuration variables is q = [\kappa_{y}], where :
    % \kappa_{y} is the curvature along y
    
    %Number of DoFs
    properties(Constant)
        n    = 1;
    end
    
    methods
        %Class constructor
        function obj = PCC2D(Parameters)
            obj      = obj@GVSBody(PCC2D.n, Parameters);
        end

        %Strain basis
        function Phi = StrainBasis(obj, s)
            Phi = [0;-1/obj.RestLength;0;0;0;0];
        end
    end

end

