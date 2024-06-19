classdef SourcePrimitive < LVPPrimitive
    %Abstract class modeling a source primitive.
    
    methods
        function obj = SourcePrimitive(Parameters)
            %Construct an instance of the stretch and compression primitive.
            arguments (Input)
                Parameters (:, 1) double
            end

            % Call the superclass constructor
            obj = obj@LVPPrimitive(Parameters);
        end
        
        % Define the methods of the primitive by overriding the default methods of the superclasses


        function [fx, dfx, ddfx, Jfq, Jfx, Jfx_ref, JJfx_q, JJfx_ref_x, JJfx_ref_q] = Update(obj, Backbone, q, dq, ddq, x, dx, ddx)
            arguments (Input)
                obj             (1, 1)  SourcePrimitive
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
                JJfx_q      (:, :)
                JJfx_ref_x  (27, :)
                JJfx_ref_q  (:, :)
            end

            %% Output preallocation for code generation
            nX              = size(x, 2);
            fx              = zeros(3, nX, "like", x);
            dfx             = zeros(3, nX, "like", x);
            ddfx            = zeros(3, nX, "like", x);
            Jfq             = zeros(3*obj.n, nX, "like", x);
            Jfx             = zeros(9, nX, "like", x);
            Jfx_ref         = zeros(9, nX, "like", x);
            JJfx_q          = zeros(9*obj.n, nX, "like", x);
            JJfx_ref_x      = zeros(27, nX, "like", x);
            JJfx_ref_q      = zeros(9*obj.n, nX, "like", x);
            
            %% Evaluate the primitive basis and its needed derivatives
            [B, dB_dx3, ddB_ddx3]             = obj.PrimitiveBasis(x(3, :));
            % Mode function
            mc                                = reshape(pagemtimes(B, q), 1, []);
            dmc_dx3                           = reshape(pagemtimes(dB_dx3, q), 1, []);
            % First order time derivative
            dx3T                              = reshape(dx(3, :), 1, 1, []);
            dB_dt                             = pagemtimes(dB_dx3, dx3T);
            dmc_dt                            = reshape(pagemtimes(B, dq) + pagemtimes(dB_dt, q), 1, []);
            % Second order time derivative
            ddx3T                             = reshape(ddx(3, :), 1, 1, []);
            ddB_ddt                           = pagemtimes(ddB_ddx3, dx3T.^2) + pagemtimes(dB_dx3, ddx3T);
            ddmc_ddt                          = reshape(pagemtimes(B, ddq) + 2*pagemtimes(dB_dt, dq) + pagemtimes(ddB_ddt, q), 1, []);

            %% Evaluate the primitive
            x1              = x(1, :);
            x2              = x(2, :);
            x3              = x(3, :);
            r_rest          = sqrt(x1.^2 + x2.^2);
            r               = sqrt(mc + x1.^2 + x2.^2);
            fx              = [(r./r_rest).*x1; (r./r_rest).*x2; x3];


            %% Evaluate the first order time derivative of the primitive
            dx1             = dx(1, :);
            dx2             = dx(2, :);
            
            x1dx1x2dx2      = x1.*dx1 + x2.*dx2;
            r_rest2         = r_rest.^2;
            a               = 1/2*(dmc_dt + 2*x1dx1x2dx2)./r;
            dfx             = [a.*x1./r_rest + r.*( dx1.*r_rest - x1.*x1dx1x2dx2./r_rest )./r_rest2; ...
                               a.*x2./r_rest + r.*( dx2.*r_rest - x2.*x1dx1x2dx2./r_rest )./r_rest2; ...
                               x(3, :)];
            
            %% Evaluate the second order time derivative of the primitive
            ddx1                = ddx(1, :);
            ddx2                = ddx(2, :);
            r2                  = r.^2;
            dmc_dt2x1dx1x2dx2   = dmc_dt + 2*x1dx1x2dx2;
            ddfx(1, :)          = (1/2)*(x1./r_rest).*( (ddmc_ddt + 2*(dx1.^2 + x1.*ddx1 + dx2.^2 + x2.*ddx2)).*r - dmc_dt2x1dx1x2dx2.*x1dx1x2dx2./r )./r2 + ...
                                + dmc_dt2x1dx1x2dx2./r.*( dx1.^2.*r_rest - x1.*x1dx1x2dx2./r_rest )./r_rest2 ...
                                + (r./(r_rest2.^2)).*( r_rest2.*( ddx1.*r_rest - (x1./r_rest2).*( dx1.^2 + x1.*ddx1 + dx2.^2 + x2.*ddx2 ).*r_rest -2*x1dx1x2dx2.^2 ) - 2*x1dx1x2dx2.*( dx1.*r_rest - x1.*x1dx1x2dx2./r_rest ) );
            ddfx(2, :)          = (1/2)*(x2./r_rest).*( (ddmc_ddt + 2*(dx1.^2 + x1.*ddx1 + dx2.^2 + x2.*ddx2)).*r - dmc_dt2x1dx1x2dx2.*x1dx1x2dx2./r )./r2 + ...
                                + dmc_dt2x1dx1x2dx2./r.*( dx2.^2.*r_rest - x2.*x1dx1x2dx2./r_rest )./r_rest2 ...
                                + (r./(r_rest2.^2)).*( r_rest2.*( ddx2.*r_rest - (x2./r_rest2).*( dx1.^2 + x1.*ddx1 + dx2.^2 + x2.*ddx2 ).*r_rest -2*x1dx1x2dx2.^2 ) - 2*x1dx1x2dx2.*( dx2.*r_rest - x2.*x1dx1x2dx2./r_rest ) );
            ddfx(3, :)          = ddx(3, :);

            %% Evaluate the Jacobian of the primitive w.r.t. x
            Jfx(1:2, :)         = [(x1.*x1)./(r.*r_rest) + r.*((r_rest - x1.*(x1./r_rest))./r_rest2); ...
                                   (x1.*x2)./(r.*r_rest) + r.*((       - x1.*(x2./r_rest))./r_rest2)];
            Jfx(4:5, :)         = [(x1.*x2)./(r.*r_rest) + r.*((       - x1.*(x2./r_rest))./r_rest2); ...
                                   (x2.*x2)./(r.*r_rest) + r.*((r_rest - x2.*(x2./r_rest))./r_rest2)];
            Jfx(7:9, :)         = [(1/2)*(dmc_dx3./r).*(x1./r_rest); ...
                                   (1/2)*(dmc_dx3./r).*(x2./r_rest); ...
                                   ones(1, nX)];
            
            
            %% Evaluate the Jacobian of the primitive w.r.t. x_ref
            % This Jacobian is zero so we do not need to do anything
            
            %% Evaluate the Jacobian of the primitive w.r.t. q
            Jfq = reshape([pagemtimes(reshape(((1/2)*(x1./r_rest)./r), 1, 1, nX), B); ...
                           pagemtimes(reshape(((1/2)*(x2./r_rest)./r), 1, 1, nX), B); ...
                           zeros(1, obj.n, nX)], 3*obj.n, nX);
            
            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x w.r.t. q
            JJfx_q = reshape((1/2).*[pagemtimes(reshape(-(x1./r_rest).*x1./(r2.^(3/2)) + (r_rest - x1.*(x1./r_rest))./(r.*r_rest2), 1, 1, []), B); ...
                                     pagemtimes(reshape(-(x2./r_rest).*x1./(r2.^(1/2)) + (       - x1.*(x2./r_rest))./(r.*r_rest2), 1, 1, []), B); ...
                                     zeros(1, obj.n, nX); ...
                                     pagemtimes(reshape(-(x1./r_rest).*x2./(r2.^(1/2)) + (       - x1.*(x2./r_rest))./(r.*r_rest2), 1, 1, []), B); ...
                                     pagemtimes(reshape(-(x2./r_rest).*x2./(r2.^(3/2)) + (r_rest - x2.*(x2./r_rest))./(r.*r_rest2), 1, 1, []), B); ...
                                     zeros(1, obj.n, nX); ...
                                     pagemtimes(reshape((x1./r_rest).*(r./r2), 1, 1, []), dB_dx3) - pagemtimes(reshape((x1./r_rest).*(dmc_dx3./r), 1, 1, []), B); ...
                                     pagemtimes(reshape((x2./r_rest).*(r./r2), 1, 1, []), dB_dx3) - pagemtimes(reshape((x2./r_rest).*(dmc_dx3./r), 1, 1, []), B); 
                                     zeros(1, obj.n, nX)], 9*obj.n, nX);

            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. x
            % This Jacobian is zero

            %% Evaluate the Jacobian of the vectorized Jacobian w.r.t. x_ref w.r.t. q
            % This Jacobian is zero
        end
    end
end

