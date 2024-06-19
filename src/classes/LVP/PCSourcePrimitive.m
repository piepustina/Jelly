classdef PCSourcePrimitive < SourcePrimitive
    %Piecewise constant source.
    
    properties
        % Number of DOFs
        n = 1
    end

    methods
        function obj = PCSourcePrimitive()
            %Construct an instance of the piecewise constant stretch and compression primitive.

            % Call the superclass
            obj = obj@SourcePrimitive([]);
        end

        % Implement the primitive basis
        function [P, dP_dx3, ddP_ddx3] = PrimitiveBasis(obj, x3)
            arguments (Input)
                obj     (1, 1) PCSourcePrimitive
                x3      (:, 1) double
            end
            arguments (Output)
                P           (1, :, :)
                dP_dx3      (1, :, :)
                ddP_ddx3    (1, :, :)
            end

            % Basis
            x3_max                  = max(x3);
            lx3                     = length(x3);
            P                       = zeros(1, obj.n, lx3);
            P(1, 1:obj.n, :)        = -(x3./x3_max).*(x3./x3_max-1);

            % Derivative of the basis w.r.t. x3
            dP_dx3                  = zeros(1, obj.n, lx3);
            dP_dx3(1, 1:obj.n, :)   = -(2*x3./(x3_max^2)-1/x3_max);

            % Second order derivative of the basis w.r.t. x3
            ddP_ddx3                = zeros(1, obj.n, lx3);
            ddP_ddx3(1, 1:obj.n, :) = -2./(x3_max^2);
        end
    end
end

