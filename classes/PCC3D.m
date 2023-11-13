classdef PCC3D < GVSBody
    %GVSBODY Class representing a slender 3D body modeled under the PCC
    %hypothesis.
    
    %Number of DoFs
    properties(Constant)
        n    = 3;
    end
    
    methods
        %Class constructor
        function obj = PCC3D(Parameters)
            obj      = obj@GVSBody(PCC3D.n, Parameters);
        end
        
        %Strain basis
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

