classdef PCC2DElongation < GVSBody
    %GVSBODY Class representing a slender body modeled under the PCC
    %hypothesis without elongation.
    
    %Number of DoFs
    properties(Constant)
        n    = 2;
    end
    
    methods
        %Class constructor
        function obj = PCC2DElongation(Parameters)
            obj      = obj@GVSBody(PCC2DElongation.n, Parameters);
        end

        %Strain basis
        function Phi = StrainBasis(obj, s)
            Phi = 1/obj.RestLength*[ 1, 0;
                                     0, 0;
                                     0, 0;
                                     0, 0;
                                     0, 0;
                                     0, 1];
        end
    end

end

