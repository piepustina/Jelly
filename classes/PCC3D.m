classdef PCC3D < GVSBody
    %Class modeling a 3D piecewise constant curvature (PCC) body with elongation.
    
    properties(Constant)
        n    = 3;
    end
    
    methods
        
        function obj = PCC3D(Parameters)
            %Construct a 3D PCC body with elongation.
            %
            %Args:
            %   Parameters ([double], [sym]): Parameters of the body, specified as for :class:`GVSBody`
            obj      = obj@GVSBody(PCC3D.n, Parameters);
        end
        
        function Phi = StrainBasis(obj, s)
            L0 = obj.RestLength;
            Phi = [   0, -1/L0,    0;
                   1/L0,     0,    0;
                      0,     0,    0;
                      0,     0,    0;
                      0,     0,    0;
                      0,     0, 1/L0];
        end
    end

end

