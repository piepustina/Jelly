classdef PAC2D < GVSBody
    %Class modeling a 2D piecewise affine curvature (PAC) body without elongation. 

    properties(Constant)
        n    = 2;
    end
    
    methods
        function obj = PAC2D(Parameters)
            %Construct a 2D PAC body without elongation.
            %
            %Args:
            %   Parameters ([double], [sym]): Parameters of the body, specified as for :class:`GVSBody`
            obj      = obj@GVSBody(PAC2D.n, Parameters);
        end

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

