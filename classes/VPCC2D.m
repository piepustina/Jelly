classdef VPCC2D < VGVSBody
    %VPCC2D Class representing a slender planar body modeled under the PCC
    %hypothesis with elongation and variable radius for two chambers. 
    
    %Number of DoFs
    properties(Constant)
        n    = 4;
    end
    
    methods
        %Class constructor
        function obj = VPCC2D(Parameters)
            obj      = obj@VGVSBody(VPCC2D.n, Parameters);
        end

        %Strain basis
        function Phi = StrainBasis(obj, s)
            Phi = 1/obj.RestLength*[ 1, 0, 0, 0;
                                     0, 0, 0, 0;
                                     0, 0, 0, 0;
                                     0, 0, 0, 0;
                                     0, 0, 0, 0;
                                     0, 1, 0, 0];
        end

        %Radius basis
        function B = RadiusBasis(obj, s, phi)
            %B = obj.GaussianRadiusBasis(s, phi);
            B = obj.BumpRadiusBasis(s, phi);
        end

        %Gaussian radius basis
        function B = GaussianRadiusBasis(obj, s, phi)
            L0  = obj.RestLength;
            if phi >= 0 && phi < pi
                %Variance
                Sigma_inv   = diag([1/(L0/2)^2; pi/2]);
                %Point evaluation
                x           = [s-L0/2; phi-pi/2];
                %Basis
                B           = [0, 0, exp(-x'*Sigma_inv*x), 0];
            else
                %Variance
                Sigma_inv   = diag([1/(L0/2)^2; pi/2]);
                %Point evaluation
                x           = [s-L0/2; phi-3*pi/2];
                %Basis
                B           = [0, 0, 0, exp(-x'*Sigma_inv*x)];
            end
        end

        %Bump radius basis
        function B = BumpRadiusBasis(obj, s, phi)
            L0  = obj.RestLength;
            if s >= L0
                B = [0, 0, 0, 0];
                return;
            end
            if phi >= 0 && phi < pi
                %Point evaluation
                x           = [s; phi];
                %Center of the bump
                x_c         = [L0/2; pi/2];
                %Basis
                B           = [0, 0, prod(exp([L0*1e-2; 1].*(-1./(x_c.^2-(x-x_c).^2)))), 0];
            else
                %Point evaluation by remapping the radius into the same
                %interval above
                x           = [s; phi-pi];
                %Center of the bump
                x_c         = [L0/2; pi/2];
                %Basis
                B           = [0, 0, 0, prod(exp([L0*1e-2; 1].*(-1./(x_c.^2-(x-x_c).^2))))];
            end
        end
        
    end

end

