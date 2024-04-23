classdef StretchCompressionPrimitive < LVPPrimitive & BackbonePrimitive
    %Abstract class modeling a stretch and compress primitive.
    
    methods
        function obj = StretchCompressionPrimitive(Parameters)
            %Construct an instance of the stretch and compression primitive.
            arguments (Input)
                Parameters (:, 1) cell
            end

            % Call the superclass constructor
            obj = obj@LVPPrimitive(Parameters);
        end
        
        % Define the methods of the primitive by overriding the default methods of the superclasses

        % Implement the strain basis function from the BackbonePrimitive superclass
        function [J, dJ] = StrainBasis(obj, s)
            arguments (Input)
                obj (1, 1) StretchCompressionPrimitive
                s   (1, 1)
            end
            
            arguments (Output)
                J   (6, :)
                dJ  (6, :)
            end

            % Output preallocation
            J   = zeros(6, obj.n);
            dJ  = zeros(6, obj.n);

            [J(6, 1:obj.n), dJ(6, 1:obj.n)] = obj.PrimitiveBasis(s);
        end

        function [fx, dfx, ddfx, Jfq, Jfx, Jfx_ref, JJfx_q, JJfx_ref_x, JJfx_ref_q] = Update(obj, Backbone, q, dq, ddq, x, dx, ddx)
            arguments (Input)
                obj             (1, 1)  StretchCompressionPrimitive
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
                Jfq         (:, :)
                Jfx         (9, :)
                Jfx_ref     (9, :)
                % Jdfx        (9, :)
                % Jdfx_ref    (9, :)
                % Jdfdx       (9, :)
                JJfx_q      (:, :)
                JJfx_ref_x  (27, :)
                JJfx_ref_q  (:, :)
            end

            %% Output preallocation for code generation
            nX              = size(x, 2);
            nBackbone       = Backbone.n;
            fx              = zeros(3, nX, "like", x);
            dfx             = zeros(3, nX, "like", x);
            ddfx            = zeros(3, nX, "like", x);
            Jfq             = zeros(3*nBackbone, nX, "like", x);
            Jfx             = zeros(9, nX, "like", x);
            Jfx_ref         = zeros(9, nX, "like", x);
            % Jdfx            = zeros(9, nX, "like", x);
            % Jdfx_ref        = zeros(9, nX, "like", x);
            % Jdfdx           = zeros(9, nX, "like", x);
            JJfx_q          = zeros(9*nBackbone, nX, "like", x);
            JJfx_ref_x      = zeros(27, nX, "like", x);
            JJfx_ref_q      = zeros(9*nBackbone, nX, "like", x);
            
            %% Evaluate the primitive
            % Get the elongation strain of the backbone, the strain is already updated by the LVPBody
            xiL               = abs(Backbone.Strain(6, 1:nX));%abs(Backbone.Strain(6, 1:nX));
            % Compute the inverse square root of the elongation strain at x
            xiL12             = xiL.^(-1/2);
            % Evaluate the primitive
            fx(1:3, 1:nX)     = [x(1:2, 1:nX).*xiL12; x(3, 1:nX) + Backbone.Elongation(1:nX)'];

            %% Evaluate the first order time derivative of the primitive
            % Extract the elongation strain
            dxiL              = Backbone.dStrain(6, 1:nX);
            % Compute the inverse square root of the elongation strain at x
            xiL32             = (xiL).^(-3/2);
            % Extract the time derivative of the elongation
            dL                = Backbone.dElongation(1:nX)';
            % Compute shared products
            xiL32dxiL = (-1/2)*xiL32.*dxiL;
            dfx(1:3, 1:nX)      = [xiL32dxiL.*x(1:2, 1:nX) + xiL12.*dx(1:2, 1:nX); ...
                                dx(3, 1:nX) + dL(1:nX)];

            %% Evaluate the second order time derivative of the primitive
            % Extract the elongation strain
            ddxiL             = Backbone.ddStrain(6, 1:nX);
            % Compute the inverse square root of the elongation strain at x
            xiL52             = (xiL).^(-5/2);
            % Extract the second order time derivative of the elongation
            ddL               = Backbone.ddElongation(1:nX)';
            % Compute commont terms
            c1                = (3/4)*xiL52.*(dxiL.^2);
            c2                = (-1/2)*xiL32.*ddxiL;
            c3                = (-1/2)*xiL32.*dxiL;
            % Compute the derivative
            ddfx(1:3, 1:nX)   = [(c1 + c2).*x(1:2, 1:nX) + 2*c3.*dx(1:2, 1:nX) + xiL12.*ddx(1:2, 1:nX); ...
                                ddx(1:nX) + ddL(1:nX)];

            %% Evaluate the Jacobian of the primitive w.r.t. x
            % NOTE: This is used in the computation of the overall Jacobian w.r.t. q and x
            Jfx(1, 1:nX)    = xiL12(1:nX);
            Jfx(5, 1:nX)    = xiL12(1:nX);
            Jfx(9, 1:nX)    = ones(1, nX, "like", x);
 
            %% Evaluate the Jacobian of the primitive w.r.t. x_ref
            xiL32Strain_ds  = (-1/2)*xiL32.*Backbone.Strain_ds(6, 1:nX);
            Jfx_ref(3, 1:nX) = xiL32Strain_ds.*x(1, 1:nX);
            Jfx_ref(6, 1:nX) = xiL32Strain_ds.*x(2, 1:nX);
            Jfx_ref(9, 1:nX) = xiL;
            
            %% Evaluate the Jacobian of the primitive w.r.t. q
            [JElongation, JElongationStrain]    = Backbone.ElongationJacobian(obj);
            JES12                               = (-1/2*xiL.^(-3/2)).*JElongationStrain;
            Jfq(1:3:3*nBackbone, 1:nX)          = JES12.*x(1, 1:nX);
            Jfq(2:3:3*nBackbone, 1:nX)          = JES12.*x(2, 1:nX);
            Jfq(3:3:3*nBackbone, 1:nX)          = JElongation;

            % %% Evaluate the Jacobian of the first order time derivative of the primitive w.r.t. x
            % Jdfx(1, 1:nX)                       = xiL32dxiL;
            % Jdfx(5, 1:nX)                       = xiL32dxiL;
            % 
            % %% Evaluate the Jacobian of the first order time derivative of the primitive w.r.t. x_ref
            % xiL_ds                              = Backbone.Strain_ds(6, 1:nX);
            % dxiL_ds                             = Backbone.dStrain_ds(6, 1:nX);
            % xiL34_52                            = (3/4)*xiL52.*dxiL.*xiL_ds + (-1/2)*xiL32.*dxiL_ds;
            % xiL32xiL_ds                         = (-1/2)*xiL32.*xiL_ds;
            % Jdfx_ref(7, 1:nX)                   = xiL34_52.*x(1, 1:nX) + xiL32xiL_ds.*dx(1, 1:nX);
            % Jdfx_ref(8, 1:nX)                   = xiL34_52.*x(2, 1:nX) + xiL32xiL_ds.*dx(2, 1:nX);
            % Jdfx_ref(9, 1:nX)                   = dxiL;
            % 
            % %% Evaluate the Jacobian of the first order time derivative of the primitive w.r.t. dx
            % Jdfdx(1, 1:nX)                      = xiL12;
            % Jdfdx(5, 1:nX)                      = xiL12;
            % Jdfdx(9, 1:nX)                      = ones(1, nX);
            
            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x w.r.t. q
            [JStrain, dJStrain]                 = Backbone.StrainJacobian(obj);
            J_xiL                               = JStrain(6, 1:nBackbone, 1:nX);
            J_xiDeltaL                          = (-1/2)*pagemtimes(reshape(xiL32, 1, 1, nX), J_xiL);
            JJfx_q                              = reshape([J_xiDeltaL; zeros(3, nBackbone, nX); J_xiDeltaL; zeros(4, nBackbone, nX)], 9*nBackbone, nX);
            
            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. x
            Zeros6x3    = zeros(6, 3, nX);
            JJfx_ref_x  = reshape([Zeros6x3; ...
                                  [reshape(xiL32Strain_ds, 1, 1, nX), zeros(1, 2, nX)]; ...
                                  [zeros(1, 1, nX) , reshape(xiL32Strain_ds, 1, 1, nX), zeros(1, 1, nX)]; ...
                                   zeros(1, 3, nX)], 27, nX);

            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. q
            Zeros                               = zeros(6, nBackbone, nX);
            JxiL_ds                             = dJStrain(6, 1:nBackbone, 1:nX);
            xiL_ds                              = Backbone.Strain_ds(6, 1:nX);
            J_xiDeltaL52x1                      = (3/4)*pagemtimes(reshape(xiL52.*xiL_ds.*x(1, 1:nX), 1, 1, nX), J_xiL);
            J_xiDeltaL_dsx1                     = (-1/2)*pagemtimes(reshape(xiL32.*x(1, 1:nX), 1, 1, nX), JxiL_ds);
            J_xiDeltaL52x2                      = (3/4)*pagemtimes(reshape(xiL52.*xiL_ds.*x(2, 1:nX), 1, 1, nX), J_xiL);
            J_xiDeltaL_dsx2                     = (-1/2)*pagemtimes(reshape(xiL32.*x(2, 1:nX), 1, 1, nX), JxiL_ds);
            JJfx_ref_q                          = reshape([Zeros; J_xiDeltaL52x1 + J_xiDeltaL_dsx1; J_xiDeltaL52x2 + J_xiDeltaL_dsx2; J_xiL], 9*nBackbone, nX);
        end

    end
end

