classdef PCC2D < GVSBody
    %Class modeling a 2D piecewise constant curvature (PCC) body without elongation.
    
    properties
        n    = 1;
    end
    
    methods
        function obj = PCC2D(Parameters)
            %Construct a 2D PCC body without elongation.
            %
            %Args:
            %   Parameters ([double], [sym]): Parameters of the body, specified as for :class:`GVSBody`
            obj      = obj@GVSBody(Parameters);
        end

        function Phi = StrainBasis(obj, s)
            Phi = [0;-1/obj.RestLength;0;0;0;0];
        end
    end

end

