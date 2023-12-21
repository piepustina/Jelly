classdef PAC3D < GVSBody
    %Class modeling a 3D piecewise affine curvature (PAC) body with elongation.

    properties(Constant)
        n    = 5;
    end
    
    methods

        function obj = PAC3D(Parameters)
            %Construct a 3D PAC body with elongation.
            %
            %Args:
            %   Parameters ([double], [sym]): Parameters of the body, specified as for :class:`GVSBody`
            obj      = obj@GVSBody(PAC3D.n, Parameters);
        end

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

