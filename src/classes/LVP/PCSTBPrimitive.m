classdef PCSTBPrimitive < STBPrimitive
    %Piecewise constant STB primitive.
    
    properties
        % Number of DOFs
        n = 5
        % Parameters for the primitive
        BodyRestLength (1, 1)
    end

    methods
        function obj = PCSTBPrimitive(BodyRestLength)
            %Construct an instance of the piecewise constant stretch and compression primitive.

            % Call the superclass
            obj = obj@STBPrimitive([BodyRestLength]);

            % Store the rest length for future use
            obj.BodyRestLength = BodyRestLength;
        end

        % Implement the primitive basis
        function [P, dP] = PrimitiveBasis(obj, x3)
            arguments (Input)
                obj (1, 1) PCSTBPrimitive
                x3  (:, 1) double
            end
            arguments (Output)
                P   (5, :, :)
                dP  (5, :, :)
            end
            

            % Output preallocation
            lx3                 = length(x3);
            P                   = zeros(5, obj.n, lx3);
            dP                  = zeros(5, obj.n, lx3);
            

            % Basis evaluation
            % Curvature
            %if x3 ~= 0 && x3 ~= obj.BodyRestLength
            %    P(1, 1)         = 1/obj.BodyRestLength;
            %    P(2, 2)         = -1/obj.BodyRestLength;
            %end
            KappaIdx                  = x3 ~= 0 & x3 ~= obj.BodyRestLength;
            P(1, 1, KappaIdx)         = 1/obj.BodyRestLength;
            P(2, 2, KappaIdx)         = -1/obj.BodyRestLength;
            
            
            % Twist
            P(3, 3, :)             = x3./obj.BodyRestLength;
            dP(3, 3, :)            = 1/obj.BodyRestLength;
            
            %Shear
            P(4, 4, :)             = x3./obj.BodyRestLength;
            P(5, 5, :)             = -x3./obj.BodyRestLength;
            dP(4, 4, :)            = 1/obj.BodyRestLength;
            dP(5, 5, :)            = -1/obj.BodyRestLength;
                        
        end
    end
end

