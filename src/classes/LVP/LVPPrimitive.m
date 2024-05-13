classdef LVPPrimitive < handle
    %LVPPRIMITIVE Abstract class that models a Locally Volume Preserving (LVP) primitive.
   
    properties(Abstract)
        %Number of DoF associated to the primitive.
        n %int32 {mustBeNonnegative}
    end

    properties (Access = public)
        % Homogeneous cell array of parameters for the primitive.
        % Each element of the parameters is a matrix of doubles that represents an argument for the contructor of the concrete classes inheriting this class.
        Parameters (:, 1) double;
    end

    properties (Constant)
        MaxParametersNumber     = 10 % Used for code generation
        MaxPrimitiveNameLength  = 50 % Used for code generation. Maximum number of characters for the primitive name
    end

    methods (Abstract)
        % Basis for the primitive
        P = PrimitiveBasis(obj, x);
    end

    methods (Static)
        function obj = loadobj(S)
            % Load a LVP Primitive
            Constructor = str2func(deblank(S.PrimitiveType));
            nParameters = sum(~isnan(S.PrimitiveParameters));
            Parameters  = S.PrimitiveParameters(1:nParameters);
            
            if isempty(Parameters)
                obj = Constructor();
            else
                obj = Constructor(Parameters);
            end
        end
    end

    methods
        function S = saveobj(obj)
            
            % Build the primitive name
            PType   = blanks(LVPPrimitive.MaxPrimitiveNameLength);
            PName   = class(obj);
            lPName  = length(PName);
            if lPName > LVPPrimitive.MaxPrimitiveNameLength
                error("For code generation, the primitive class name cannot exceed " + LVPPrimitive.MaxPrimitiveNameLength + " characters.");
            end 
            PType(1:lPName) = PName;

            % Build the primitive parameters
            PParameters  = nan(1, LVPPrimitive.MaxParametersNumber);
            lP = length(obj.Parameters);
            if lP > LVPPrimitive.MaxParametersNumber
                error("For code generation, the number of primitive parameters cannot exceed " + LVPPrimitive.MaxParametersNumber + " .");
            end
            PParameters(1:lP) = obj.Parameters';

            S = struct('PrimitiveType', PType, 'PrimitiveParameters', PParameters);
        end
    end
    
    methods
        function obj = LVPPrimitive(Parameters)
            %Construct the locally volume preserving primitve.

            arguments (Input)
                Parameters (:, 1) double
            end

            arguments (Output)
                obj (1, 1) LVPPrimitive
            end
            
            % Store the parameters for the primitive
            obj.Parameters = Parameters;
        end

        function S = toStruct(obj)
            % Convert the object primitive to a struct
            S = struct('PrimitiveType', class(obj), 'PrimitiveParameters', {obj.Parameters});
        end
        
        % Update method for the primitive
        function [fx, dfx, ddfx, Jfq, Jfx, Jfx_ref, JJfx_q, JJfx_ref_x, JJfx_ref_q] = Update(obj, Backbone, q, dq, ddq, x, dx, ddx)
            arguments (Input)
                obj             (1, 1)  LVPPrimitive
                Backbone        (1, 1)  LVPBackbone
                q               (:, 1)                  = zeros(obj.n, 1)
                dq              (:, 1)                  = zeros(obj.n, 1)
                ddq             (:, 1)                  = zeros(obj.n, 1)
                x               (3, :)                  = zeros(3, 1)
                dx              (3, :)                  = zeros(3, 1)
                ddx             (3, :)                  = zeros(3, 1)
            end

            arguments (Output)
                fx          (3, :)
                dfx         (3, :)
                ddfx        (3, :)
                % Jacobian of the primitive w.r.t. q
                Jfq         (:, :)
                % Jacobian of the primitive w.r.t. x.
                % NOTE: This Jacobian is used to compute both the Jacobian w.r.t. q and x (underformed)
                Jfx         (9, :)
                % Jacobian of the primitive w.r.t. x_ref where x_ref is the reference value of x
                % NOTE: This Jacobian is used, together with Jfxq, to compute the Jacobian w.r.t. x.
                Jfx_ref     (9, :)
                % % Jacobian of the time derivative of the primitive w.r.t. x
                % Jdfx        (9, :)
                % % Jacobian of the time derivative of the primitive w.r.t. x_ref
                % Jdfx_ref    (9, :)
                % % Jacobian of the time derivative of the primitive w.r.t dx
                % Jdfdx       (9, :)
                % Jacobian of the vectorized Jacobian w.r.t. x w.r.t. q
                JJfx_q      (:, :)
                % Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. x
                JJfx_ref_x  (27, :)
                % Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. q
                JJfx_ref_q  (:, :)
            end

            % By default a primitive does not transform the input
            Nx              = size(x, 2);
            nBackbone       = Backbone.n;
            fx              = x;
            dfx             = dx;
            ddfx            = ddx;
            Jfq             = zeros(3*nBackbone, Nx, "like", x);
            Jfx             = zeros(9, Nx, "like", x);
            Jfx_ref         = zeros(9, Nx, "like", x);
            % Jdfx            = zeros(9, Nx, "like", x);
            % Jdfx_ref        = zeros(9, Nx, "like", x);
            % Jdfdx           = zeros(9, Nx, "like", x);
            JJfx_q          = zeros(9*nBackbone, Nx, "like", x);
            JJfx_ref_x      = zeros(27, Nx, "like", x);
            JJfx_ref_q      = zeros(9*nBackbone, Nx, "like", x);
            Ones            = ones(1, Nx, "like", x);
            Jfx(1, 1:Nx)    = Ones;
            Jfx(5, 1:Nx)    = Ones;
            Jfx(9, 1:Nx)    = Ones;
        end
        
    end
end

