classdef Dual
    %DUAL Define a dual vector for automatic differentiation.
    
    properties
        x  (:, :) double
        dx (:, :) double
        Nx (1, 1) int32
    end

    methods (Static)
        % Check if x is a dual or a double
        function bool = isDualOrDouble(x)
            bool = Dual.isDual(x) || isa(x, "double");
        end

        function bool = isDual(x)
            bool = isa(x, "Dual");
        end

        % Cast from double to dual
        function d = DoubleToDual(x)
            % Check if the variable is already a dual
            if Dual.isDual(x)
                d = Dual(x.x, x.dx);
            else
                d = Dual(x);
            end
        end
    end
    
    methods
        function obj = Dual(x, isVariable)
            %DUAL Construct a dual object
            arguments
                x          (:, :) {Dual.isDualOrDouble}
                isVariable (1, 1) = true
            end
            % Check if the first argument is a dual object or not
            if isa(x, "Dual")
                obj = x;
            else
                % Check that x is a vector
                if ~isvector(x)
                    error("x must be a vector");
                end
                obj.Nx = length(x);
                obj.x  = x;
                % Assign the derivative
                if isVariable == true
                    obj.dx= eye(length(x));
                else
                    obj.dx= zeros(length(x));
                end
            end
        end

        % Custom operator for indexing
        function c = slice(a, idx)
            arguments (Input)
                a (1, 1) {Dual.isDual}
                idx (:, :) int32
            end
            arguments (Output)
                c (1, 1) Dual
            end
            c    = Dual(a.x(idx));
            c.dx = a.dx(idx, :);
        end
        
        %% Overload MATLAB operators
        % Sum
        function c = plus(a, b)
            arguments (Input)
                a (:, :) {Dual.isDualOrDouble}
                b (:, :) {Dual.isDualOrDouble}
            end
            arguments (Output)
                c (1, 1) Dual
            end
            % Cast a and b to dual
            ad = Dual(a, false);
            bd = Dual(b, false);
            % Compute the output
            c    = Dual(ad.x + bd.x);
            c.dx = ad.dx + bd.dx;
        end
        
        % Difference
        function c = minus(a, b)
            arguments (Input)
                a (:, :) {Dual.isDualOrDouble}
                b (:, :) {Dual.isDualOrDouble}
            end
            arguments (Output)
                c (1, 1) Dual
            end
            % Cast a and b to dual
            ad = Dual(a, false);
            bd = Dual(b, false);
            % Compute the output
            c    = Dual(ad.x - bd.x);
            c.dx = ad.dx - bd.dx;
        end
        
        % Unitary sum
        function c = uplus(a)
            arguments (Input)
                a (:, :) {Dual.isDualOrDouble}
            end
            arguments (Output)
                c (1, 1) Dual
            end
            c = Dual(a, false);
        end
        
        % Unitary difference
        function c = uminus(a)
            arguments (Input)
                a (:, :) {Dual.isDualOrDouble}
            end
            arguments (Output)
                c (1, 1) Dual
            end
            % Cast a to Dual
            ad = Dual(a, false);
            % Compute the output
            c = Dual(-ad.x);
            c.dx = -ad.dx;
        end
        
        % Element wise product
        function c = times(a, b)
            arguments (Input)
                a (:, :) {Dual.isDualOrDouble}
                b (:, :) {Dual.isDualOrDouble}
            end
            arguments (Output)
                c (1, 1) Dual
            end
            % Cast a and b to dual
            ad = Dual(a, false);
            bd = Dual(b, false);
            % Compute the output
            c = Dual(ad.x .* bd.x);
            c.dx = ad.dx*bd.x + ad.x*bd.dx;
        end
        
        % Matrix product
        function c = mtimes(a, b)
            arguments (Input)
                a (:, :) {Dual.isDualOrDouble}
                b (:, :) {Dual.isDualOrDouble}
            end
            arguments (Output)
                c (1, 1) Dual
            end
            % Cast a and b to dual
            ad = Dual(a, false);
            bd = Dual(b, false);
            % Compute the output
            c = Dual(ad.x * bd.x);
            c.dx = ad.dx*bd.x + ad.x*bd.dx;
        end
        
        % Element-wise right division
        function c = rdivide(a, b)
            arguments (Input)
                a (:, :) {Dual.isDualOrDouble}
                b (:, :) {Dual.isDualOrDouble}
            end
            arguments (Output)
                c (1, 1) Dual
            end
            % Cast a and b to dual
            ad = Dual(a, false);
            bd = Dual(b, false);
            % Compute the output
            c = Dual(ad.x / bd.x, (ad.dx*bd.x - ad.x*bd.dx)/bd.x^2);
        end
        
        % Element-wise power
        function c = power(a, b)
            arguments
                a (:, :) {Dual.isDualOrDouble}
                b (:, :) {Dual.isDualOrDouble}
            end
            % Cast a and b to dual
            ad = Dual(a);
            bd = Dual(b);
            c = Dual(ad.x.^bd.x, bd.x*ad.x.^(bd.x-1)*ad.dx);
        end

        % Matrix power
        function c = mpower(a, b)
            arguments
                a (:, :) {Dual.isDualOrDouble}
                b (1, 1) {Dual.isDualOrDouble}
            end
            % Cast a and b to dual
            ad = Dual(a);
            bd = Dual(b);
            c = Dual(ad.x^bd.x, bd.x*ad.x^(bd.x-1)*ad.dx);
        end

        % Matrix transpose
        function at = transpose(a)
            arguments (Input)
                a (1, 1) {Dual.isDual}
            end
            arguments (Output)
                at (1, 1) Dual
            end
            at   = Dual(a);
            at.x = a.x.';
        end
        
        % Complex conjugate tranpose
        function at = ctranspose(a)
            arguments (Input)
                a (1, 1) {Dual.isDual}
            end
            arguments (Output)
                at (1, 1) Dual
            end
            at   = Dual(a);
            at.x = a.x';
        end
        
        % Display function
        function disp(a)
            disp("Value: ");
            disp(a.x);
            disp("Derivative: ");
            disp(a.dx);
        end

    end
end

