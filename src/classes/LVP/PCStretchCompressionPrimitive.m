classdef PCStretchCompressionPrimitive < StretchCompressionPrimitive
    %Piecewise constant stretch and compression primitive.
    
    properties
        % Number of DOFs
        n = 1
        % Parameters for the primitive
        BodyRestLength (1, 1)
    end

    methods
        function obj = PCStretchCompressionPrimitive(BodyRestLength)
            %Construct an instance of the piecewise constant stretch and compression primitive.

            % Call the superclass
            obj = obj@StretchCompressionPrimitive({BodyRestLength});

            % Store the rest length for future use
            obj.BodyRestLength = BodyRestLength;
        end

        % Implement the primitive basis
        function [P, dP] = PrimitiveBasis(obj, x3)
            arguments (Input)
                obj (1, 1) PCStretchCompressionPrimitive
                x3  (1, 1) double
            end
            arguments (Output)
                P   (1, :)
                dP  (1, :)
            end

            % Basis
            P           = zeros(1, obj.n);
            P(1:obj.n)  = (x3/obj.BodyRestLength)*(x3/obj.BodyRestLength-1);

            % Derivative of the basis w.r.t. x3
            dP          = zeros(1, obj.n);
            dP(1:obj.n) = 2*x3/(obj.BodyRestLength^2)-1/obj.BodyRestLength;
        end
    end
end

