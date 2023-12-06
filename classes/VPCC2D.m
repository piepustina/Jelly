classdef VPCC2D < VGVSBody
    %VPCC2D Class representing a slender planar body modeled under the PCC
    %hypothesis with elongation and variable radius.
    
    %Number of DoFs
    properties(Constant)
        n    = 3;
    end
    
    methods
        %Class constructor
        function obj = VPCC2D(Parameters)
            obj      = obj@VGVSBody(VPCC2D.n, Parameters);
        end

        %Strain basis
        function Phi = StrainBasis(obj, s)
            Phi = 1/obj.RestLength*[ 0, 0, 0;
                                    -1, 0, 0;
                                     0, 0, 0;
                                     0, 0, 0;
                                     0, 0, 0;
                                     0, 1, 0];
        end

        %Radius basis
        function B = RadiusBasis(obj, s)
            B = [0, 0, 1];
        end
    end

end

