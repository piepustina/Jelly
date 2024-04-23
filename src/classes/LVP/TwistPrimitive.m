classdef TwistPrimitive < LVPPrimitive & BackbonePrimitive
    %Class modeling a twist primitive
    
    
    methods
        % Class constructor
        function obj = TwistPrimitive(Parameters)
            %Construct an instance of the twist primitive.
            arguments (Input)
                Parameters (:, 1) cell
            end

            % Call the superclass constructor
            obj = obj@LVPPrimitive(Parameters);
        end

        % Implement the strain basis function from the BackbonePrimitive superclass
        function [J, dJ] = StrainBasis(obj, s)
            arguments (Input)
                obj (1, 1) TwistPrimitive
                s   (1, 1)
            end

            arguments (Output)
                J   (6, :)
                dJ  (6, :)
            end

            % Output preallocation
            J   = zeros(6, obj.n);
            dJ  = zeros(6, obj.n);
            
            % Call the primitive basis
            [J(3, 1:obj.n), dJ(3, 1:obj.n)] = obj.PrimitiveBasis(s);
        end

        % Update method for the primitive
        function [fx, dfx, ddfx, Jfq, Jfx, Jfx_ref, JJfx_q, JJfx_ref_x, JJfx_ref_q] = Update(obj, Backbone, q, dq, ddq, x, dx, ddx)
            arguments (Input)
                obj             (1, 1)  TwistPrimitive
                Backbone        (1, 1)  LVPBackbone
                q               (:, 1)                  = zeros(obj.n, 1)
                dq              (:, 1)                  = zeros(obj.n, 1)
                ddq             (:, 1)                  = zeros(obj.n, 1)
                x               (3, :)                  = zeros(3, 1)
                dx              (3, :)                  = zeros(3, 1)
                ddx             (3, :)                  = zeros(3, 1)
            end

            arguments (Output)
                fx              (3, :)
                dfx             (3, :)
                ddfx            (3, :)
                Jfq             (:, :)
                Jfx             (9, :)
                Jfx_ref         (9, :)
                % Jdfx            (9, :)
                % Jdfx_ref        (9, :)
                % Jdfdx           (9, :)
                JJfx_q          (:, :)
                JJfx_ref_x      (27, :)
                JJfx_ref_q      (:, :)
            end

            %% Output preallocation
            Nx              = size(x, 2);
            nBackbone       = Backbone.n;
            fx              = zeros(3, Nx, "like", x);
            dfx             = zeros(3, Nx, "like", x);
            ddfx            = zeros(3, Nx, "like", x);
            Jfq             = zeros(3*nBackbone, Nx, "like", x);
            Jfx             = zeros(9, Nx, "like", x);
            Jfx_ref         = zeros(9, Nx, "like", x);
            % Jdfx            = zeros(9, Nx, "like", x);
            % Jdfx_ref        = zeros(9, Nx, "like", x);
            % Jdfdx           = zeros(9, Nx, "like", x);
            JJfx_q          = zeros(9*nBackbone, Nx, "like", x);
            JJfx_ref_x      = zeros(27, Nx, "like", x);
            JJfx_ref_q      = zeros(9*nBackbone, Nx, "like", x);

            %% Evaluate the primitive
            % Get the pose of the backbone at the query points 
            Pose              = Backbone.Pose;
            % Compute the position in the global frame of each backbone points
            t                 = squeeze(Pose(1:3, 4, 1:Nx));
            n1                = squeeze(Pose(1:3, 1, 1:Nx));
            n2                = squeeze(Pose(1:3, 2, 1:Nx));
            % Evaluate the primitive
            fx(1:3, 1:Nx)     = t + x(1, 1:Nx).*n1 + x(2, 1:Nx).*n2;

            %% Evaluate the first order time derivative
            omega             = Backbone.dPose(1:3, 1:Nx);
            dt                = Backbone.dPose(4:6, 1:Nx);
            dn1               = cross(omega, n1);
            dn2               = cross(omega, n2);
            dfx(1:3, 1:Nx)    = dt + ...
                                dx(1, 1:Nx).*n1 + x(1, 1:Nx).*dn1 + ...
                                dx(2, 1:Nx).*n2 + x(2, 1:Nx).*dn2;

            %% Evaluate the second order time derivative
            domega            = Backbone.ddPose(1:3, 1:Nx);
            ddt               = Backbone.ddPose(4:6, 1:Nx);
            ddn1              = cross(domega, n1) + cross(omega, dn1);
            ddn2              = cross(domega, n2) + cross(omega, dn2);
            ddfx(1:3, 1:Nx)   = ddt + ...
                                ddx(1, 1:Nx).*n1 + 2*dx(1, 1:Nx).*dn1 + x(1, 1:Nx).*ddn1 + ...
                                ddx(2, 1:Nx).*n2 + 2*dx(2, 1:Nx).*dn2 + x(2, 1:Nx).*ddn2;
            
            %% Evaluate the Jacobian w.r.t. x
            Jfx(1:3, 1:Nx)    = n1(1:3, 1:Nx);
            Jfx(4:6, 1:Nx)    = n2(1:3, 1:Nx);

            %% Evaluate the Jacobian w.r.t. x_ref
            AngularStrain       = Backbone.Strain(1:3, 1:Nx);
            LinearStrain        = Backbone.Strain(4:6, 1:Nx);
            R                   = Pose(1:3, 1:3, 1:Nx);
            dR_dx3              = pagemtimes(R, skew(AngularStrain));% Derivative of R (rotation matrix) w.r.t. x3_ref
            dn1_dx3             = squeeze(dR_dx3(1:3, 1, 1:Nx));
            dn2_dx3             = squeeze(dR_dx3(1:3, 2, 1:Nx));
            dt_dx3              = squeeze(pagemtimes(R, reshape(LinearStrain, 3, 1, Nx)));% Derivative of the origin position w.r.t. x3_ref
            Jfx_ref(7:9, 1:Nx)  = dt_dx3 + x(1, 1:Nx).*dn1_dx3 + x(2, 1:Nx).*dn2_dx3;

            %% Evaluate the Jacobian w.r.t. q
            JPose             = Backbone.PoseJacobian(obj);
            JL                = reshape(JPose(4:6, 1:nBackbone, 1:Nx), 3*nBackbone, Nx);
            JA                = JPose(1:3, 1:nBackbone, 1:Nx);
            Jn1               = reshape(pagemtimes(-1*skew(n1), JA), 3*nBackbone, Nx);
            Jn2               = reshape(pagemtimes(-1*skew(n2), JA), 3*nBackbone, Nx);
            Jfq               = JL + x(1, 1:Nx).*Jn1 + x(2, 1:Nx).*Jn2;

            % %% Evaluate the Jacobian of the first order time derivative of the primitive w.r.t. x
            % Jdfx(1:3, 1:Nx)     = dn1;
            % Jdfx(4:6, 1:Nx)     = dn2;
            % 
            % %% Evaluate the Jacobian of the first order time derivative of the primitive w.r.t. x_ref
            % dAngularStrain      = Backbone.dStrain(1:3, 1:Nx);
            % dLinearStrain       = Backbone.dStrain(4:6, 1:Nx);
            % 
            % AngularStrainG      = reshape(pagemtimes(R, reshape(AngularStrain, 3, 1, Nx)), 3, Nx);
            % LinearStrainG       = pagemtimes(R, reshape(LinearStrain, 3, 1, Nx));
            % dAngularStrainG     = reshape(pagemtimes(R, reshape(dAngularStrain, 3, 1, Nx)), 3, Nx);
            % dLinearStrainG      = pagemtimes(R, reshape(dLinearStrain, 3, 1, Nx));
            % 
            % dn1_dt_ds           = pagemtimes(skew(omega), pagemtimes(skew(AngularStrainG), Pose(1:3, 1, 1:Nx))) + pagemtimes(skew(dAngularStrainG), Pose(1:3, 1, 1:Nx)); 
            % dn2_dt_ds           = pagemtimes(skew(omega), pagemtimes(skew(AngularStrainG), Pose(1:3, 2, 1:Nx))) + pagemtimes(skew(dAngularStrainG), Pose(1:3, 2, 1:Nx)); 
            % dt_dt_ds            = pagemtimes(skew(omega), LinearStrainG) + dLinearStrainG;
            % 
            % Jdfx_ref(7:9, 1:Nx) = squeeze(dt_dt_ds) ...
            %                     + dn1_dx3.*dx(1, 1:Nx) + squeeze(dn1_dt_ds).*x(1, 1:Nx)...
            %                     + dn2_dx3.*dx(2, 1:Nx) + squeeze(dn2_dt_ds).*x(2, 1:Nx);
            % 
            % %% Evaluate the Jacobian of the first order time derivative of the primitive w.r.t. dx
            % Jdfdx(1:3, 1:Nx)    = n1;
            % Jdfdx(4:6, 1:Nx)    = n2;

            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x w.r.t. q
            Zeros  = zeros(3, nBackbone, Nx);
            JJfx_q = reshape([pagemtimes(-1*skew(n1), JA); pagemtimes(-1*skew(n2), JA); Zeros], 9*nBackbone, Nx);
            
            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. x
            Zeros6x3    = zeros(6, 3, Nx);
            JJfx_ref_x  = reshape([Zeros6x3; reshape(dn1_dx3, 3, 1, Nx), reshape(dn2_dx3, 3, 1, Nx), zeros(3, 1, Nx)], 27, Nx);

            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. q
            Jxi                                 = Backbone.StrainJacobian(obj);
            % Angular contributions
            JxiA                                = Jxi(1:3, 1:nBackbone, 1:Nx);
            SkewR_xiA                           = skew(reshape(-1*pagemtimes(R, reshape(AngularStrain, 3, 1, Nx)), 3, Nx));
            R_JxiA                              = -1*pagemtimes(R, JxiA);
            Jn1_q                               = pagemtimes(skew(reshape(pagemtimes(SkewR_xiA, Pose(1:3, 1, 1:Nx)), 3, Nx)), JA) + pagemtimes(skew(n1), R_JxiA);
            Jn2_q                               = pagemtimes(skew(reshape(pagemtimes(SkewR_xiA, Pose(1:3, 2, 1:Nx)), 3, Nx)), JA) + pagemtimes(skew(n2), R_JxiA);
            % Linear contributions
            SkewR_xiL                           = skew(reshape(-1*pagemtimes(R, reshape(LinearStrain, 3, 1, Nx)), 3, Nx));
            JxiL                                = Jxi(4:6, 1:nBackbone, 1:Nx);
            Jt_q                                = pagemtimes(SkewR_xiL, JPose(1:3, 1:nBackbone, 1:Nx)) + pagemtimes(R, JxiL);
            % Overall Jacobian
            JJfx_ref_q                          = reshape([Zeros; Zeros; Jt_q + reshape(x(1, 1:Nx), 1, 1, Nx).*Jn1_q + reshape(x(2, 1:Nx), 1, 1, Nx).*Jn2_q], 9*nBackbone, Nx);
        end
    end
end

