classdef PCTwistShearPrimitive < TwistShearPrimitive
    %Piecewise constant shear and twist primitive.
    
    properties
        % Number of DOFs
        n = 3
        % Parameters for the primitive
        BodyRestLength (1, 1)
    end

    methods
        function obj = PCTwistShearPrimitive(BodyRestLength)
            %Construct an instance of the piecewise constant stretch and compression primitive.

            % Call the superclass
            obj = obj@TwistShearPrimitive([BodyRestLength]);

            % Store the rest length for future use
            obj.BodyRestLength = BodyRestLength;
        end

        % Implement the primitive basis
        function [P, dP] = PrimitiveBasis(obj, x3)
            arguments (Input)
                obj (1, 1) PCTwistShearPrimitive
                x3  (1, 1) double
            end
            arguments (Output)
                P  (3, :)
                dP (3, :)
            end
            
            % Basis 
            P                 = zeros(3, obj.n);
            P(1:3, 1:obj.n)   = [x3/obj.BodyRestLength,                     0,                     0;
                                                     0, x3/obj.BodyRestLength,                     0; 
                                                     0,                     0,-x3/obj.BodyRestLength];

            % Derivative of the basis w.r.t. x3 
            dP                = zeros(3, obj.n);
            dP(1:3, 1:obj.n)  = [1/obj.BodyRestLength,                    0,                    0;
                                                    0, 1/obj.BodyRestLength,                    0; 
                                                    0,                    0,-1/obj.BodyRestLength];
        end
    end
end

