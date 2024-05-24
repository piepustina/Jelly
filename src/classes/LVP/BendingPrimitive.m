classdef BendingPrimitive < LVPPrimitive & BackbonePrimitive 
    %Class that models a bending primitive
    
    methods
        % Class constructor
        function obj = BendingPrimitive(Parameters)
            %Construct an instance of the bending primitive.
            arguments (Input)
                Parameters (:, 1) double
            end

            % Call the superclass constructor
            obj = obj@LVPPrimitive(Parameters);
        end

        % Implement the strain basis function and its derivative w.r.t. s from the BackbonePrimitive superclass
        function [J, dJ] = StrainBasis(obj, s)
            arguments (Input)
                obj (1, 1) BendingPrimitive
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
            [J(1:2, 1:obj.n), dJ(1:2, 1:obj.n)] = obj.PrimitiveBasis(s);
        end

        % Update method for the primitive
        function [fx, dfx, ddfx, Jfq, Jfx, Jfx_ref, JJfx_q, JJfx_ref_x, JJfx_ref_q] = Update(obj, Backbone, q, dq, ddq, x, dx, ddx)
            arguments (Input)
                obj             (1, 1)  BendingPrimitive
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
            % Cartesian position of x
            x1                = x(1, 1:Nx);
            x2                = x(2, 1:Nx);
            

            % Get the curvature strains
            Kappa_x           = Backbone.Strain(1, 1:Nx);
            Kappa_y           = Backbone.Strain(2, 1:Nx);
            % Get the elongation strain
            Elongation        = Backbone.Strain(6, 1:Nx);
            
            % Compute the polar coordinates of the points x1 and x2
            [theta, Rad]      = cart2pol(x1, x2);
            % Evaluate sin and cos
            c_theta           = cos(theta);
            s_theta           = sin(theta);
            % Compute the coefficients for the radius restriction
            a                 = (2/3)*(Kappa_x.*s_theta - Kappa_y.*c_theta)./(Elongation);
            a                 = abs(a);% We do consider the absolute value
            c                 = Rad;
            % Find the index of points for which a ~= 0
            idxAisZero        = a <= 1e-4; % Threshold for the computation of a and to run limits
            % Preallocate the radius expansion
            a_c               = complex(a);% Use complex values for the computation of the radius and for code generation
            c_c               = complex(c);
            a2c2              = (a_c.^2).*(c_c.^2);
            d                 = (27*a2c2 + sqrt((27*a2c2-2).^2 - 4) - 2).^(1/3);
            r                 = (real(d./(3*2^(1/3)*a_c) + (2^(1/3))./(3*a_c.*d) - 1./(3*a_c)));
            r(idxAisZero)     = Rad(idxAisZero);
            % Compute the primitive
            fx(1:3, 1:Nx)     = t + r.*c_theta.*n1 + r.*s_theta.*n2;
            % Get the index of the points with zero radius
            idx               = Rad == 0;

            %% Evaluate the first order time derivative of the primitive
            % Partial derivatives of theta w.r.t. x1 and x2
            dx1_dt            = dx(1, 1:Nx);
            dx2_dt            = dx(2, 1:Nx);
            dtheta_dx1        = -x2./(x1.^2 + x2.^2);
            dtheta_dx2        =  x1./(x1.^2 + x2.^2);
            dtheta_dx1(idx)   = 0;
            dtheta_dx2(idx)   = 0;
            % First order time derivative of theta
            dtheta_dt         = dtheta_dx1.*dx1_dt + dtheta_dx2.*dx2_dt;
            
            % Partial derivatives of r w.r.t. a and c parameters
            dr_da             = -(r.^2)./(2 + 3*a.*r);
            dr_dc             = 2*c./(r.*(2 + 3*a.*r));
            dr_dc(idx)        = 0;
            % First order time derivative of the curvature strains
            dKappa_x_dt       = Backbone.dStrain(1, 1:Nx);
            dKappa_y_dt       = Backbone.dStrain(2, 1:Nx);
            % First order time derivative of the elongation
            dElongation_dt    = Backbone.dStrain(6, 1:Nx);

            % First order time derivative of a and c
            %da_dt             = (2/3)*( dKappa_x_dt.*s_theta + Kappa_x.*c_theta.*dtheta_dt - dKappa_y_dt.*c_theta + Kappa_y.*s_theta.*dtheta_dt );
            da_dt             = (2/3)*1./(Elongation.^2).*(Elongation.*( dKappa_x_dt.*s_theta + Kappa_x.*c_theta.*dtheta_dt - dKappa_y_dt.*c_theta + Kappa_y.*s_theta.*dtheta_dt ) - (Kappa_x.*s_theta - Kappa_y.*c_theta).*dElongation_dt);
            
            dc_dt             = (x1.*dx1_dt + x2.*dx2_dt)./sqrt(x1.^2 + x2.^2);
            dc_dt(idx)        = 0;
            % First order time derivative of r
            dr_dt             = dr_da.*da_dt + dr_dc.*dc_dt;
            % Angular velocitity of the backbone
            omega             = Backbone.dPose(1:3, 1:Nx);
            % First order time derivative of position vector t
            dt_dt             = Backbone.dPose(4:6, 1:Nx);
            % First order time derivative of n1 and n2
            dn1_dt            = cross(omega, n1);
            dn2_dt            = cross(omega, n2);
            % Compute the time derivative of the primitive
            dfx(1:3, 1:Nx)    = dt_dt ...
                              + (dr_dt.*c_theta - r.*s_theta.*dtheta_dt).*n1 + r.*c_theta.*dn1_dt ...
                              + (dr_dt.*s_theta - r.*c_theta.*dtheta_dt).*n2 + r.*s_theta.*dn2_dt;

            %% Evaluate the second order time derivative
            % Compute the second order time derivative of theta
            ddx1_ddt          = ddx(1, 1:Nx);
            ddx2_ddt          = ddx(2, 1:Nx);
            ddtheta_ddt       = ( (x1.^2 + x2.^2).*( x1.*ddx2_ddt - x2.*ddx1_ddt ) + (x1.*dx2_dt - x2.*dx1_dt).*2.*(x1 .* dx1_dt + x2.*dx2_dt) )./(( x1.^2 + x2.^2 ).^2);
            ddtheta_ddt(idx)  = 0;
            % Compute the first order time derivative of the partial derivative of dr_da and dr_dc
            ddr_da_dt         = -( 2*r.*dr_dt.*( 2 + 3*a.*r ) - (r.^2).*( 3*da_dt.*r + 3*a.*dr_dt ) )./((2 + 3*a.*r).^2);
            ddr_dc_dt         = ( 2*dc_dt.*( r.*( 2 + 3*a.*r ) ) - 2*c.*( dr_dt.*( 2 + 3*a.*r ) + 3*r.*( da_dt.*r + a.*dr_dt ) ) )./( (r.^2).*( 2 + 3*a.*r ).^2 );
            ddr_dc_dt(idx)    = 0;
            % Second order time derivative of the curvature strains
            ddKappa_x_ddt     = Backbone.ddStrain(1, 1:Nx);
            ddKappa_y_ddt     = Backbone.ddStrain(2, 1:Nx);
            % Second order time derivative of the elongation strain
            ddElongation_ddt  = Backbone.ddStrain(6, 1:Nx);

            % Compute the second order time derivative of a and c
            % dda_ddt           = (2/3)*(( ddKappa_x_ddt.*s_theta + 2*dKappa_x_dt.*c_theta.*dtheta_dt - Kappa_x.*s_theta.*dtheta_dt.^2 + Kappa_x.*c_theta.*ddtheta_ddt ...
            %                             -ddKappa_y_ddt.*c_theta + 2*dKappa_y_dt.*s_theta.*dtheta_dt + Kappa_y.*c_theta.*dtheta_dt.^2 + Kappa_y.*s_theta.*ddtheta_ddt));
            % 
            dda_ddt           = (2/3)*1./(Elongation.^4).*( ((Elongation.^2).*(Elongation.*( ddKappa_x_ddt.*s_theta + 2*dKappa_x_dt.*c_theta.*dtheta_dt - Kappa_x.*s_theta.*dtheta_dt.^2 + Kappa_x.*c_theta.*ddtheta_ddt ...
                                                            -ddKappa_y_ddt.*c_theta + 2*dKappa_y_dt.*s_theta.*dtheta_dt + Kappa_y.*c_theta.*dtheta_dt.^2 + Kappa_y.*s_theta.*ddtheta_ddt) - (Kappa_x.*s_theta - Kappa_y.*c_theta)).*ddElongation_ddt) - ( Elongation.*(dKappa_x_dt.*s_theta + Kappa_x.*c_theta.*dtheta_dt - dKappa_y_dt.*c_theta + Kappa_y.*s_theta.*dtheta_dt)- (Kappa_x.*s_theta - Kappa_y.*c_theta).*dElongation_dt).*(2*Elongation).*dElongation_dt );
            
            ddc_ddt           = ( sqrt(x1.^2 + x2.^2).*( dx1_dt.^2 + x1.*ddx1_ddt + dx2_dt.^2 + x2.*ddx2_ddt) - (x1.*dx1_dt + x2.*dx2_dt).*(x1.*dx1_dt + x2.*dx2_dt)./(sqrt( x1.^2 + x2.^2 )) )./( x1.^2 + x2.^2 );
            ddc_ddt(idx)      = 0;

            % Compute the second order time derivative of the r
            ddr_ddt           = ddr_da_dt.*da_dt + dr_da.*dda_ddt + ddr_dc_dt.*dc_dt + dr_dc.*ddc_ddt;
            
            % Angular acceleration
            domega_dt         = Backbone.ddPose(1:3, 1:Nx);
            % Second order time derivative of the backbone position vector
            ddt_ddt           = Backbone.ddPose(4:6, 1:Nx);
            % Second order time derivative of n1 and n2
            ddn1_dt           = cross(domega_dt, n1) + cross(omega, dn1_dt);
            ddn2_dt           = cross(domega_dt, n2) + cross(omega, dn2_dt);
            
            % Second order time derivative of the primitive
            ddfx(1:3, 1:Nx)   = ddt_ddt ...
                              + ddr_ddt.*c_theta.*n1 - dr_dt.*s_theta.*dtheta_dt.*n1 + dr_dt.*c_theta.*dn1_dt ...
                              - dr_dt.*s_theta.*dtheta_dt.*n1 - r.*c_theta.*dtheta_dt.^2.*n1 - r.*s_theta.*ddtheta_ddt.*n1 - r.*s_theta.*dtheta_dt.*dn1_dt ...
                              + dr_dt.*c_theta.*dn1_dt - r.*s_theta.*dtheta_dt.*dn1_dt + r.*c_theta.*ddn1_dt ...
                              + ddr_ddt.*s_theta.*n2 + dr_dt.*c_theta.*dtheta_dt.*n2 + dr_dt.*s_theta.*dn2_dt ...
                              + dr_dt.*c_theta.*dtheta_dt.*n2 - r.*s_theta.*dtheta_dt.^2.*n2 + r.*c_theta.*ddtheta_ddt.*n2 + r.*c_theta.*dtheta_dt.*dn2_dt ...
                              + dr_dt.*s_theta.*dn2_dt + r.*c_theta.*dtheta_dt.*dn2_dt + r.*s_theta.*ddn2_dt;
            
            %% Evaluate the Jacobian w.r.t. x
            % Derivative of a and c w.r.t. x1
            %da_dx1              = (2/3)*( Kappa_x.*c_theta.*dtheta_dx1 + Kappa_y.*s_theta.*dtheta_dx1 );
            da_dx1              = (2/3)*(1./Elongation).*( Kappa_x.*c_theta.*dtheta_dx1 + Kappa_y.*s_theta.*dtheta_dx1 );
            
            dc_dx1              = x1./sqrt( x1.^2 + x2.^2 );
            dc_dx1(idx)         = 0;
            % Derivative of a and c w.r.t. x2
            %da_dx2              = (2/3)*( Kappa_x.*c_theta.*dtheta_dx2 + Kappa_y.*s_theta.*dtheta_dx2 );
            da_dx2              = (2/3)*(1./Elongation).*( Kappa_x.*c_theta.*dtheta_dx2 + Kappa_y.*s_theta.*dtheta_dx2 );
            
            dc_dx2              = x2./sqrt( x1.^2 + x2.^2 );
            dc_dx2(idx)         = 0;
            % Derivative of r w.r.t. x1
            dr_dx1              = dr_da.*da_dx1 + dr_dc.*dc_dx1;
            % Derivative of r w.r.t. x2
            dr_dx2              = dr_da.*da_dx2 + dr_dc.*dc_dx2;
            
            % Jacobian w.r.t x1
            Jfx(1:3, 1:Nx)      = dr_dx1.*c_theta.*n1 - r.*s_theta.*dtheta_dx1.*n1 + dr_dx1.*s_theta.*n2 + r.*c_theta.*dtheta_dx1.*n2;
            % Jacobian w.r.t. x2
            Jfx(4:6, 1:Nx)      = dr_dx2.*c_theta.*n1 - r.*s_theta.*dtheta_dx2.*n1 + dr_dx2.*s_theta.*n2 + r.*c_theta.*dtheta_dx2.*n2;

            %% Evaluate the Jacobian w.r.t. x_ref
            % Derivative of the curvature strains w.r.t. x3_ref
            dKappa_x_dx3_ref    = Backbone.dStrain_ds(1, 1:Nx);
            dKappa_y_dx3_ref    = Backbone.dStrain_ds(2, 1:Nx);
            % Derivative of the elongation strain w.r.t. x3_ref
            dElongation_dx3_ref = Backbone.dStrain_ds(6, 1:Nx);
            
            % Derivative of a w.r.t. x3_ref
            %da_dx3_ref          = (2/3)*( dKappa_x_dx3_ref.*s_theta - dKappa_y_dx3_ref.*c_theta );
            da_dx3_ref          = (2/3)*(1./(Elongation.^2)).*(Elongation.*( dKappa_x_dx3_ref.*s_theta - dKappa_y_dx3_ref.*c_theta ) - (Kappa_x.*s_theta - Kappa_y.*c_theta).*dElongation_dx3_ref);
            
            % Derivative of r w.r.t. x3_ref
            dr_dx3_ref          = dr_da.*da_dx3_ref;
            
            % Compute the derivatives w.r.t x3_ref of the backbone quantities
            AngularStrain       = Backbone.Strain(1:3, 1:Nx);
            LinearStrain        = Backbone.Strain(4:6, 1:Nx);
            R                   = Pose(1:3, 1:3, 1:Nx);
            % Derivative of the rotation matrix w.r.t. x3_ref
            dR_dx3_ref          = pagemtimes(R, skew(AngularStrain));% Derivative of R (rotation matrix) w.r.t. x3_ref
            % Derivative of n1 and n2 w.r.t. x3_ref
            dn1_dx3_ref         = squeeze(dR_dx3_ref(1:3, 1, 1:Nx));
            dn2_dx3_ref         = squeeze(dR_dx3_ref(1:3, 2, 1:Nx));
            % Derivative of the position vector t w.r.t. x3_ref
            dt_dx3_ref          = squeeze(pagemtimes(R, reshape(LinearStrain, 3, 1, Nx)));% Derivative of the origin position w.r.t. x3_ref
            % Jacobian w.r.t. x3_ref
            Jfx_ref(7:9, 1:Nx)  = dt_dx3_ref + dr_dx3_ref.*c_theta.*n1 + r.*c_theta.*dn1_dx3_ref + dr_dx3_ref.*s_theta.*n2 + r.*s_theta.*dn2_dx3_ref;

            %% Evaluate the Jacobian w.r.t. q
            % Get the Jacobian of the strain w.r.t. q and that of the derivative w.r.t. x3_ref w.r.t. q 
            [JStrain_dq, JdStrain_dx3_refdq] = Backbone.StrainJacobian(obj);
            % Derivative of the curvatures w.r.t. q
            dKappa_x_dq       = squeeze(JStrain_dq(1, 1:nBackbone, 1:Nx));
            dKappa_y_dq       = squeeze(JStrain_dq(2, 1:nBackbone, 1:Nx));
            % Derivatives of the elongation w.r.t. q
            dElongation_dq    = squeeze(JStrain_dq(6, 1:nBackbone, 1:Nx));

            % Derivative of a w.r.t. q
            %da_dq             = (2/3)*( dKappa_x_dq.*s_theta - dKappa_y_dq.*c_theta );
            da_dq             = (2/3)*1./(Elongation.^2).*(Elongation.*( dKappa_x_dq.*s_theta - dKappa_y_dq.*c_theta ) - (Kappa_x.*s_theta - Kappa_y.*c_theta).*dElongation_dq);
            % Derivative of r w.r.t. q
            dr_dq             = dr_da.*da_dq;

            % Derivatives of the backbone w.r.t. q
            JPose             = Backbone.PoseJacobian(obj);
            Jt_dq             = reshape(JPose(4:6, 1:nBackbone, 1:Nx), 3*nBackbone, Nx);
            JA                = JPose(1:3, 1:nBackbone, 1:Nx);
            Jn1_dq            = reshape(pagemtimes(-skew(n1), JA), 3*nBackbone, Nx);
            Jn2_dq            = reshape(pagemtimes(-skew(n2), JA), 3*nBackbone, Nx);
            
            % Evaluate the primitive Jacobian w.r.t. q
            Jfq               = Jt_dq ...
                              + reshape(pagemtimes(reshape(dr_dq, 1, 1, nBackbone, Nx), repmat(reshape(c_theta.*n1, 3, 1, 1, Nx), 1, 1, nBackbone, 1)), 3*nBackbone, Nx) + r.*c_theta.*Jn1_dq ...
                              + reshape(pagemtimes(reshape(dr_dq, 1, 1, nBackbone, Nx), repmat(reshape(s_theta.*n2, 3, 1, 1, Nx), 1, 1, nBackbone, 1)), 3*nBackbone, Nx) + r.*s_theta.*Jn2_dq;

            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x w.r.t. q
            % Derivative of dr_dc w.r.t. q
            ddr_dcdq          = (-2*c./(r.^2.*( 2 + 3*a.*r).^2)).*( dr_dq.*( 2 + 3*a.*r) + 3*r.*( da_dq.*r + a.*dr_dq ) );

            % Derivative of da_dx1 w.r.t.q 
            %dda_dx1dq         = (2/3)*( dKappa_x_dq.*c_theta.*dtheta_dx1 + dKappa_y_dq.*s_theta.*dtheta_dx1 );
            dda_dx1dq         = (2/3)*(1./(Elongation.^2)).*( Elongation.*( dKappa_x_dq.*c_theta.*dtheta_dx1 + dKappa_y_dq.*s_theta.*dtheta_dx1 ) - (Kappa_x.*c_theta.*dtheta_dx1 + Kappa_y.*s_theta.*dtheta_dx1).*dElongation_dq );

            % Derivative of da_dx2 w.r.t.q 
            %dda_dx2dq         = (2/3)*( dKappa_x_dq.*c_theta.*dtheta_dx2 + dKappa_y_dq.*s_theta.*dtheta_dx2 );
            dda_dx2dq         = (2/3)*(1./(Elongation.^2)).*( Elongation.*( dKappa_x_dq.*c_theta.*dtheta_dx2 + dKappa_y_dq.*s_theta.*dtheta_dx2 ) - (Kappa_x.*c_theta.*dtheta_dx2 + Kappa_y.*s_theta.*dtheta_dx2).*dElongation_dq );

            % Derivative of dr_da w.r.t. q
            ddr_dadq          = -( 2*r.*dr_dq.*( 2 + 3*a.*r ) - 3*(r.^2).*( da_dq.*r + a.*dr_dq ) )./( (2 + 3*a.*r).^2 );

            % Derivative of dr_dx1 w.r.t. q
            ddr_dx1dq         = ddr_dadq.*da_dx1 + dr_da.*dda_dx1dq + ddr_dcdq.*dc_dx1;

            % Derivative of dr_dx2 w.r.t. q
            ddr_dx2dq         = ddr_dadq.*da_dx2 + dr_da.*dda_dx2dq + ddr_dcdq.*dc_dx2;

            % Derivative of n1 and n2 w.r.t.q
            dn1_dq            = pagemtimes(-skew(n1), JA);
            dn2_dq            = pagemtimes(-skew(n2), JA);

            % Rows 1-3 of JJfx_q
            JJfx1_q= squeeze(pagemtimes(reshape(ddr_dx1dq, 1, 1, nBackbone, Nx), repmat(reshape(c_theta.*n1 + s_theta.*n2, 3, 1, 1, Nx), 1, 1, nBackbone, 1))) ...
                   + squeeze(pagemtimes(reshape(dr_dq, 1, 1, nBackbone, Nx)    , repmat(reshape(-s_theta.*dtheta_dx1.*n1 + c_theta.*dtheta_dx1.*n2, 3, 1, 1, Nx), 1, 1, nBackbone, 1))) ...
                   + pagemtimes(reshape(dr_dx1.*c_theta - r.*s_theta.*dtheta_dx1, 1, 1, Nx), dn1_dq) ...
                   + pagemtimes(reshape(dr_dx1.*s_theta + r.*c_theta.*dtheta_dx1, 1, 1, Nx), dn2_dq);

            % Rows 4-6 of JJfx_q
            JJfx2_q= squeeze(pagemtimes(reshape(ddr_dx2dq, 1, 1, nBackbone, Nx), repmat(reshape(c_theta.*n1 + s_theta.*n2, 3, 1, 1, Nx), 1, 1, nBackbone, 1))) ...
                   + squeeze(pagemtimes(reshape(dr_dq, 1, 1, nBackbone, Nx)    , repmat(reshape(-s_theta.*dtheta_dx2.*n1 + c_theta.*dtheta_dx2.*n2, 3, 1, 1, Nx), 1, 1, nBackbone, 1))) ...
                   + pagemtimes(reshape(dr_dx2.*c_theta - r.*s_theta.*dtheta_dx2, 1, 1, Nx), dn1_dq) ...
                   + pagemtimes(reshape(dr_dx2.*s_theta + r.*c_theta.*dtheta_dx2, 1, 1, Nx), dn2_dq);

            Zeros  = zeros(3, nBackbone, Nx);
            JJfx_q = reshape([JJfx1_q; JJfx2_q; Zeros], 9*nBackbone, Nx);
            
            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. x
            % Derivative of dr_da w.r.t. x1 and x2
            ddr_dadx1   = -( 2*r.*dr_dx1.*( 2 +3*a.*r ) - (3*r.^2).*( da_dx1.*r + a.*dr_dx1 ) )./((2 + 3*a.*r).^2);
            ddr_dadx2   = -( 2*r.*dr_dx2.*( 2 +3*a.*r ) - (3*r.^2).*( da_dx2.*r + a.*dr_dx2 ) )./((2 + 3*a.*r).^2);

            % Derivative of da_dx3_ref w.r.t. x1 and x2
            %dda_dx3_refdx1 = (2/3)*( dKappa_x_dx3_ref.*c_theta.*dtheta_dx1 + dKappa_y_dx3_ref.*s_theta.*dtheta_dx1 );
            %dda_dx3_refdx2 = (2/3)*( dKappa_x_dx3_ref.*c_theta.*dtheta_dx2 + dKappa_y_dx3_ref.*s_theta.*dtheta_dx2 );
            dda_dx3_refdx1 = (2/3)*(1./(Elongation.^2)).*(Elongation.*( dKappa_x_dx3_ref.*c_theta.*dtheta_dx1 + dKappa_y_dx3_ref.*s_theta.*dtheta_dx1 ) - dElongation_dx3_ref.*(Kappa_x.*c_theta.*dtheta_dx1 + Kappa_y.*s_theta.*dtheta_dx1) );
            dda_dx3_refdx2 = (2/3)*(1./(Elongation.^2)).*(Elongation.*( dKappa_x_dx3_ref.*c_theta.*dtheta_dx2 + dKappa_y_dx3_ref.*s_theta.*dtheta_dx2 ) - dElongation_dx3_ref.*(Kappa_x.*c_theta.*dtheta_dx2 + Kappa_y.*s_theta.*dtheta_dx2) );


            % Derivative of dr_dx3_ref w.r.t. x1 and x2
            ddr_dx3_refdx1 = ddr_dadx1.*da_dx3_ref + dr_da.*dda_dx3_refdx1;
            ddr_dx3_refdx2 = ddr_dadx2.*da_dx3_ref + dr_da.*dda_dx3_refdx2;
                       
            % Compute the first column of JJfx_ref_x
            JJfx_ref_x3_dx1 = ddr_dx3_refdx1.*c_theta.*n1 - dr_dx3_ref.*s_theta.*dtheta_dx1.*n1 + dr_dx1.*c_theta.*dn1_dx3_ref - r.*s_theta.*dtheta_dx1.*dn1_dx3_ref ...
                            + ddr_dx3_refdx1.*s_theta.*n2 + dr_dx3_ref.*c_theta.*dtheta_dx1.*n2 + dr_dx1.*s_theta.*dn2_dx3_ref + r.*c_theta.*dtheta_dx1.*dn2_dx3_ref;
            JJfx_ref_x3_dx2 = ddr_dx3_refdx2.*c_theta.*n1 - dr_dx3_ref.*s_theta.*dtheta_dx2.*n1 + dr_dx2.*c_theta.*dn1_dx3_ref - r.*s_theta.*dtheta_dx2.*dn1_dx3_ref ...
                            + ddr_dx3_refdx2.*s_theta.*n2 + dr_dx3_ref.*c_theta.*dtheta_dx2.*n2 + dr_dx2.*s_theta.*dn2_dx3_ref + r.*c_theta.*dtheta_dx2.*dn2_dx3_ref;

            Zeros6x3    = zeros(6, 3, Nx);
            JJfx_ref_x  = reshape([Zeros6x3; reshape(JJfx_ref_x3_dx1, 3, 1, Nx), reshape(JJfx_ref_x3_dx2, 3, 1, Nx), zeros(3, 1, Nx)], 27, Nx);

            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. q

            % Derivative of dKappa_x_dx3_ref and dKappa_y_dx3_ref w.r.t. q
            ddKappa_x_dx3_refdq = squeeze(JdStrain_dx3_refdq(1, 1:nBackbone, 1:Nx));
            ddKappa_y_dx3_refdq = squeeze(JdStrain_dx3_refdq(2, 1:nBackbone, 1:Nx));

            % Derivative of dElongation_dx3_ref w.r.t. q
            ddElongation_dx3_refdq = squeeze(JdStrain_dx3_refdq(6, 1:nBackbone, 1:Nx));

            % Derivative of da_dx3_ref w.r.t. q
            %dda_dx3_refdq       = (2/3)*( ddKappa_x_dx3_refdq.*s_theta - ddKappa_y_dx3_refdq.*c_theta );
            dda_dx3_refdq       = (2/3)*(1./(Elongation.^4)).*( (Elongation.^2).*( dElongation_dq.*(dKappa_x_dx3_ref.*s_theta-dKappa_y_dx3_ref.*c_theta) + Elongation.*(ddKappa_x_dx3_refdq.*s_theta - ddKappa_y_dx3_refdq.*c_theta) - ( dElongation_dx3_ref.*( dKappa_x_dq.*s_theta - dKappa_y_dq.*c_theta ) + ( Kappa_x.*s_theta - Kappa_y.*c_theta ).*ddElongation_dx3_refdq )) - 2*Elongation.*dElongation_dq.*( Elongation.*(dKappa_x_dx3_ref.*s_theta - dKappa_y_dx3_ref.*c_theta) - (Kappa_x.*s_theta - Kappa_y.*c_theta).*dElongation_dx3_ref ) );
            
            % Derivative of dr_dx3_ref w.r.t. q
            ddr_dx3_refdq       = ddr_dadq.*da_dx3_ref + dr_da.*dda_dx3_refdq;

            % Jacobian of dn1_dx3_ref and dn2_dx3_ref w.r.t. q
            JxiA                                = JStrain_dq(1:3, 1:nBackbone, 1:Nx);
            SkewR_xiA                           = skew(reshape(pagemtimes(R, reshape(AngularStrain, 3, 1, Nx)), 3, Nx));
            Skewn1                              = skew(n1);
            Skewn2                              = skew(n2);
            R_JxiA                              = pagemtimes(R, JxiA);
            Jdn1_dx3ref_dq                      = pagemtimes(pagemtimes(Skewn1, SkewR_xiA)-pagemtimes(SkewR_xiA, Skewn1), JA) - pagemtimes(Skewn1, R_JxiA);
            Jdn2_dx3ref_dq                      = pagemtimes(pagemtimes(Skewn2, SkewR_xiA)-pagemtimes(SkewR_xiA, Skewn2), JA) - pagemtimes(Skewn2, R_JxiA);
            
            % Jacobian of dt_dx3_ref w.r.t. q
            SkewR_xiL                           = skew(reshape(pagemtimes(R, reshape(LinearStrain, 3, 1, Nx)), 3, Nx));
            JxiL                                = JStrain_dq(4:6, 1:nBackbone, 1:Nx);
            Jdt_dx3_ref_dq                      = -pagemtimes(SkewR_xiL, JA) + pagemtimes(R, JxiL);
            
            % Last three rows of JJfx_ref_q
            JJfx_ref_q3                         = Jdt_dx3_ref_dq ...
                                                + squeeze(pagemtimes(reshape(ddr_dx3_refdq, 1, 1, nBackbone, Nx), repmat(reshape(c_theta.*n1, 3, 1, 1, Nx), 1, 1, nBackbone, 1)))...
                                                + pagemtimes(reshape(dr_dx3_ref.*c_theta, 1, 1, Nx), dn1_dq) ...
                                                + squeeze(pagemtimes(reshape(dr_dq, 1, 1, nBackbone, Nx), repmat(reshape(c_theta.*dn1_dx3_ref, 3, 1, 1, Nx), 1, 1, nBackbone, 1)))...
                                                + pagemtimes(reshape(r.*c_theta, 1, 1, Nx), Jdn1_dx3ref_dq) ...
                                                + squeeze(pagemtimes(reshape(ddr_dx3_refdq, 1, 1, nBackbone, Nx), repmat(reshape(s_theta.*n2, 3, 1, 1, Nx), 1, 1, nBackbone, 1)))...
                                                + pagemtimes(reshape(dr_dx3_ref.*s_theta, 1, 1, Nx), dn2_dq) ...
                                                + squeeze(pagemtimes(reshape(dr_dq, 1, 1, nBackbone, Nx), repmat(reshape(s_theta.*dn2_dx3_ref, 3, 1, 1, Nx), 1, 1, nBackbone, 1)))...
                                                + pagemtimes(reshape(r.*s_theta, 1, 1, Nx), Jdn2_dx3ref_dq);

            % Overall Jacobian
            JJfx_ref_q                          = reshape([Zeros; Zeros; JJfx_ref_q3], 9*nBackbone, Nx);
        end
    end

    methods (Access = protected)
        % Compute the radius
        function R = Radius(obj, a_r, c_r)
            arguments (Input)
                obj (1, 1) BendingPrimitive
                a_r   (1, :) double 
                c_r   (1, :) double
            end
            arguments (Output)
                R   (1, :) double
            end
            % Preallocate the output
            R_c = complex(c_r);
            a   = complex(a_r);
            c   = complex(c_r);
            
            % Compute the new radius if a is larger than the limit threshold
            idxa = a>= 1e-5;
            R_c(1, idxa) = (-1./(2.*a(1, idxa)) + ...
                    3./(2.^(2/3).*a(1, idxa).*(-54 + 324*(a(1, idxa).^2).*(c(1, idxa).^2) + sqrt(-2916 + (-54 + 324*(a(1, idxa).^2).*(c(1, idxa).^2)).^2)).^(1/3)) +...
                    (-54 + 324*(a(1, idxa).^2).*(c(1, idxa).^2) + sqrt(-2916 + (-54 + 324*(a(1, idxa).^2).*(c(1, idxa).^2)).^2)).^(1/3)./(6*2^(1/3).*a(1, idxa)));
            % Take the real part as the above operation operates on complex
            R = real(R_c);
        end
    end
end

