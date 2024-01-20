classdef PCC2DElongation < GVSBody
    %Class modeling a 2D piecewise constant curvature (PCC) body with elongation.

    properties(Constant)
        n    = 2;
    end
    
    methods
        
        function obj = PCC2DElongation(Parameters)
            %Construct a 2D PCC body with elongation.
            %
            %Args:
            %   Parameters ([double], [sym]): Parameters of the body, specified as for :class:`GVSBody`
            obj      = obj@GVSBody(PCC2DElongation.n, Parameters);
        end

        function Phi = StrainBasis(obj, s)
            Phi = 1/obj.RestLength*[ 0, 0;
                                    -1, 0;
                                     0, 0;
                                     0, 0;
                                     0, 0;
                                     0, 1];
        end
    end

end

