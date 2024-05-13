classdef PCTwistPrimitive < TwistPrimitive
    %Piecewise constant twist primitive.
    
    properties
        % Number of DOFs
        n = 1
        % Parameters for the primitive
        BodyRestLength (1, 1)
    end

    methods
        function obj = PCTwistPrimitive(BodyRestLength)
            %Construct an instance of the piecewise constant stretch and compression primitive.

            % Call the superclass
            obj = obj@TwistPrimitive([BodyRestLength]);

            % Store the rest length for future use
            obj.BodyRestLength = BodyRestLength;
        end

        % Implement the primitive basis
        function [P, dP] = PrimitiveBasis(obj, x3)
            arguments (Input)
                obj (1, 1) PCTwistPrimitive
                x3  (1, 1) double
            end
            arguments (Output)
                P   (1, :)
                dP  (1, :)
            end
            
            % Basis evaluation
            P               = zeros(1, obj.n);
            P(1:obj.n)      = x3/obj.BodyRestLength;
            
            % Dervivative of the basis w.r.t. x3
            dP              = zeros(1, obj.n);
            dP(1:obj.n)     = 1/obj.BodyRestLength;
        end
    end
end

