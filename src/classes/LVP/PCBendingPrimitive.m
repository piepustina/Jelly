classdef PCBendingPrimitive < BendingPrimitive
    %Piecewise constant bending primitive.
    
    properties
        % Number of DOFs
        n = 2
        % Parameters for the primitive
        BodyRestLength (1, 1)
    end

    methods
        function obj = PCBendingPrimitive(BodyRestLength)
            %Construct an instance of the piecewise constant stretch and compression primitive.

            % Call the superclass
            obj = obj@BendingPrimitive([BodyRestLength]);

            % Store the rest length for future use
            obj.BodyRestLength = BodyRestLength;
        end

        % Implement the primitive basis
        function [P, dP] = PrimitiveBasis(obj, x3)
            arguments (Input)
                obj (1, 1) PCBendingPrimitive
                x3  (1, 1) double
            end
            arguments (Output)
                P   (2, :)
                dP  (2, :)
            end
            
            % Basis evaluation
            P                   = zeros(2, obj.n);

            P(1:2, 1:obj.n)     = [1/obj.BodyRestLength,  0; ...
                                                      0, -1/obj.BodyRestLength];

            % P(1:2, 1:obj.n)     = [sin(pi*x3/obj.BodyRestLength),  0; ...
            %                                           0, -sin(pi*x3/obj.BodyRestLength)];

            % Derivative of the basis w.r.t. x3
            dP                  = zeros(2, obj.n);
            % dP(1:2, 1:obj.n)    = [cos(pi*x3/obj.BodyRestLength)*pi/obj.BodyRestLength,  0; ...
            %                                           0, -cos(pi*x3/obj.BodyRestLength)*pi/obj.BodyRestLength];
        end
    end
end

