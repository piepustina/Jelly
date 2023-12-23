classdef GVSBody < Body
    %Abstract class modeling a 1D continuum under the Geometric Variable Strain (GVS) approach. The internal interaction forces are modeled using a linear visco-elastic law.
    
    methods (Abstract)
        %Strain basis as a function of the arc length
        StrainBasis(obj, s);
    end

    properties
        %Vector collecting the parameters of the GVS Body as :math:`\mathrm{Parameters} = \left( L_{0} \,\, R_{\mathrm{base}} \,\, R_{\mathrm{tip}} \,\, \rho \,\, E \,\, \nu \,\, \eta \,\, N_{\mathrm{Gauss}} \right)^{T} \in \mathbb{R}^{8 \times 1}`, where :math:`L_{0}, R_{\mathrm{base}}` and :math:`R_{\mathrm{tip}}` denote the body rest length, the base and tip radius, :math:`\rho` is the mass density, :math:`E` the Young modulus, :math:`\nu`
        %the Poisson ratio, :math:`\eta` the material damping coefficient, and :math:`N_{\mathrm{Gauss}}` represents the number of Gaussian points used for the computation of the kinematics and the integrals.
        Parameters;
        RestLength      = 0;
        ReferenceStrain = [0;0;0;0;0;1];
    end

    % Private properties extracted from the parameters for ease of use
    properties (Access = private)
        BaseRadius      = 0;
        TipRadius       = 0;
        MassDensity     = 0;
        YoungModulus    = 0;
        ShearModulus    = 0;
        PoissonRatio    = 0.5;
        MaterialDamping = 0;
        NGaussPoints    = 0;
        GaussPoints     = 0;
        GaussWeights    = 0;
        NGaussPointsInt = 0;
        % Properties used for the internal computation, all the quantities are expressed in the global frame
        gGauss;
        EtaGauss;
        JEtaGauss;
        dEtaGauss;
        rGauss;
        drGauss;
        JdrGauss;
        ddrGauss;
    end
    
    methods
        %% Class constructor
        function obj = GVSBody(n, Parameters)
            %Construct a GVS body.
            %
            %Args:
            %   n          (double)         : Number of DoF of the body
            %   Parameters ([double], [sym]): Parameters of the body, specified as :math:`L_{0}, R_{\mathrm{base}}, R_{\mathrm{tip}}, \rho, E, \nu, \eta` and :math:`N_{\mathrm{Gauss}}`
            
            %Convert the parameters to a column vector, if needed
            if isrow(Parameters)
                Parameters = Parameters';
            end
            % Call the superclass constructor
            obj                     = obj@Body(n);
            % Set the parameters
            obj.Parameters          = Parameters;
            obj.RestLength          = Parameters(1);
            obj.BaseRadius          = Parameters(2);
            obj.TipRadius           = Parameters(3);
            obj.MassDensity         = Parameters(4);
            obj.YoungModulus        = Parameters(5);
            obj.PoissonRatio        = Parameters(6);
            obj.ShearModulus        = obj.YoungModulus/(2*(1+obj.PoissonRatio));
            obj.MaterialDamping     = Parameters(7);
            obj.NGaussPoints        = Parameters(8);
            % Preallocate the Gaussian points
            [obj.GaussPoints, obj.GaussWeights] = lgwt(obj.NGaussPoints, 0, obj.RestLength);
            obj.GaussPoints         = [obj.GaussPoints; obj.RestLength];
            obj.GaussWeights        = [obj.GaussWeights; 0];
            obj.NGaussPointsInt     = obj.NGaussPoints + 1;
            obj.gGauss              = repmat(eye(4)     , [1, 1, obj.NGaussPointsInt]);
            obj.EtaGauss            = zeros(6, obj.NGaussPointsInt);
            obj.JEtaGauss           = zeros(6, obj.n, obj.NGaussPointsInt);
            obj.dEtaGauss           = zeros(6, obj.NGaussPointsInt);
            obj.rGauss              = zeros(3, obj.NGaussPointsInt);
            obj.drGauss             = zeros(3, obj.NGaussPointsInt);
            obj.JdrGauss            = zeros(3, obj.n, obj.NGaussPointsInt);
            obj.ddrGauss            = zeros(3, obj.NGaussPointsInt);
        end

        % Radius function
        function R = Radius(obj, s)
            R = obj.TipRadius + (obj.RestLength - s)/obj.RestLength*(obj.BaseRadius - obj.TipRadius);
        end

        % Mass density along the curvliniear abscissa
        function rhos = LinearMassDensity(obj, s)
            rhos = pi*obj.Radius(s)^2*obj.MassDensity;
        end

        %% Body methods implementation
        % Overload the update body method
        function Update(obj, q, dq, ddq)
            [T_, omega_, v_, domega_, dv_, J_omega, J_v] = obj.Kinematics(q, dq, ddq); % Update the kinematics of the body
            obj.T_                  = T_;
            obj.v_rel_              = v_;
            obj.omega_rel_          = omega_;
            obj.a_rel_              = dv_;
            obj.domega_rel_         = domega_;
            RT                      = T_(1:3, 1:3)';
            %The jacobians are expressed in the tip frame
            obj.v_par_              = RT*J_v;
            obj.omega_par_          = RT*J_omega;
            %Inertial quantities, the configuration variables are not
            %required since all the computations are done using internal
            %state variables
            obj.p_com_              = obj.p_com();
            obj.v_com_rel_          = obj.v_com_rel();
            obj.a_com_rel_          = obj.a_com_rel();
            obj.I_                  = obj.I();
            obj.m_                  = obj.m();
            obj.J_                  = obj.J();
            obj.int_dr_             = obj.int_dr(q, dq);%TODO: To be removed since it is not required
            obj.int_ddr_            = obj.int_ddr(q, dq, ddq);%TODO: To be removed since it is always zero
            obj.int_r_X_dr_         = obj.int_r_X_dr();
            obj.int_r_X_ddr_        = obj.int_r_X_ddr();
            obj.int_dr_X_pv_r_      = obj.int_dr_X_pv_r();
            obj.int_pv_r_O_dd_r_    = obj.int_pv_r_O_dd_r();
            obj.int_dr_O_dr_        = obj.int_dr_O_dr(q, dq);%TODO: To be removed since it is always zero
            obj.grad_int_dr_        = obj.grad_int_dr(q);%TODO: To be removed since it is always zero
            obj.grad_int_r_X_dr_    = obj.grad_int_r_X_dr();
            obj.grad_J_             = obj.grad_J();
            obj.grad_v_com_         = obj.grad_v_com();
            obj.K_                  = obj.K(q);
            obj.D_                  = obj.D(q, dq);
        end
        % Transformation from base to body tip
        function T_ = T(obj, q)
            T_ = obj.T_s(q, obj.RestLength);
        end
        % Transformation matrix from base to s
        function Ts_ = T_s(obj, q, s)
            [Ts_, ~, ~, ~, ~, ~, ~] = obj.Kinematics_s(q, zeros(obj.n, 1, 'like', q), zeros(obj.n, 1, 'like', q), s);
        end
        % Relative velocity of the body tip
        function v_rel_ = v_rel(obj, q, dq)
            [~, ~, v_rel_, ~, ~, ~, ~] = obj.Kinematics(q, dq, zeros(obj.n, 1, 'like', q));
        end
        % Relative angular velocity
        function omega_rel_ = omega_rel(obj, q, dq)
            [~, omega_rel_, ~, ~, ~, ~, ~] = obj.Kinematics(q, dq, zeros(obj.n, 1, 'like', q));
        end
        % Relative linear acceleration of the tip
        function a_rel_ = a_rel(obj, q, dq, ddq)
            [~, ~, ~, ~, a_rel_, ~, ~] = obj.Kinematics(q, dq, ddq);
        end
        % Relative angular acceleration
        function domega_rel_ = domega_rel(obj, q, dq, ddq)
            [~, ~, ~, domega_rel_, ~, ~, ~] = obj.Kinematics(q, dq, ddq);
        end
        % Jacobian of the linear velocity of the tip with respect to q in the tip frame
        function v_par_ = v_par(obj, q)
            [g, ~, ~, ~, ~, ~, v_par_] = obj.Kinematics(q, zeros(obj.n, 1, 'like', q), zeros(obj.n, 1, 'like', q));
            v_par_ = g(1:3, 1:3)'*v_par_;
        end
        % Jacobian of the angular velocity with respect to q in the tip frame
        function omega_par_ = omega_par(obj, q)
            [g, ~, ~, ~, ~, omega_par_, ~] = obj.Kinematics(q, zeros(obj.n, 1, 'like', q), zeros(obj.n, 1, 'like', q));
            omega_par_ = g(1:3, 1:3)'*omega_par_;
        end
        % Center of mass position
        function p_com_  = p_com(obj, q)
            if nargin == 2 % Update the DK
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end
            % Return the CoM position comuted in Kinematics
            p_com_  = obj.p_com_;
        end
        % Relative velocity of the center of mass
        function v_com_rel_ = v_com_rel(obj, q, dq)
            if nargin == 2 % Update the DK
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            else
                if nargin == 3
                    obj.Kinematics(q, dq, zeros([obj.n, 1], 'like', q));
                end
            end
            % Return the CoM velocity
            v_com_rel_   = obj.v_com_rel_;
        end
        % Relative acceleration of the center of mass
        function a_com_rel_  = a_com_rel(obj, q, dq, ddq)
            if nargin == 2 % Update the DK
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            else
                if nargin == 3
                    obj.Kinematics(q, dq, zeros([obj.n, 1], 'like', q));
                else
                    if nargin == 4
                        obj.Kinematics(q, dq, ddq);
                    end
                end
            end
            % Compute the acceleration
            a_com_rel_ = obj.a_com_rel_;
        end
        
        % Inertia matrix
        function I_ = I(obj, q)
            if nargin == 2 % Update the DK
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end
            % Compute the inertia
            I_ = zeros(3, 3);
            for i = 1:obj.NGaussPointsInt
                Sr_i = skew(obj.rGauss(:, i));
                I_  = I_ + obj.GaussWeights(i)*(Sr_i'*Sr_i*obj.LinearMassDensity(obj.GaussPoints(i)));
            end
        end
        
        % Mass
        function m_ = m(obj)
            m_ = 0;
            for i = 1:obj.NGaussPointsInt
                m_ = m_ + obj.LinearMassDensity(obj.GaussPoints(i))*obj.GaussWeights(i);
            end
        end
        
        % Time derivative of the inertia matrix
        function J_ = J(obj, q, dq)
            if nargin == 2%Update the kinematics
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            else
                if nargin == 3
                    obj.Kinematics(q, dq, zeros([obj.n, 1], 'like', q));
                end
            end
            %Compute the time derivative of the inertia
            J_ = zeros(3, 3, 'like', obj.rGauss(1:3, 1));
            for i = 1:obj.NGaussPointsInt
                Sr_i  = skew(obj.rGauss(:, i));
                Sdr_i = skew(obj.drGauss(:, i));
                J_ = J_ + (Sdr_i'*Sr_i + Sr_i'*Sdr_i)*obj.LinearMassDensity(obj.GaussPoints(i))*obj.GaussWeights(i);
            end
        end
        
        % Integral of \dot{r}
        % TODO: This is always zero, remove
        function int_dr_ = int_dr(obj, q, dq)
            int_dr_ = zeros(3, 1, 'like', q);
        end

        % Integral of \ddot{r}
        % TODO: This is always zero, remove
        function int_ddr_ = int_ddr(obj, q, dq, ddq)
            int_ddr_ = zeros(3, 1, 'like', q);
        end

        % Integral of \cross(r, \dot{r})
        function int_r_X_dr_ = int_r_X_dr(obj, q, dq)
            if nargin == 2% Update the kinematics
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            else
                if nargin == 3
                    obj.Kinematics(q, dq, zeros([obj.n, 1], 'like', q));
                end
            end
            % Compute the term
            int_r_X_dr_ = zeros(3, 1, 'like', obj.rGauss(1:3, 1));
            for i = 1:obj.NGaussPointsInt
                int_r_X_dr_ = int_r_X_dr_ + cross(obj.rGauss(:, i), obj.drGauss(:, i))*obj.LinearMassDensity(obj.GaussPoints(i))*obj.GaussWeights(i);
            end
        end

        % Integral of \cross(r, \ddot{r})
        function int_r_X_ddr_ = int_r_X_ddr(obj, q, dq, ddq)
            if nargin == 2% Update the kinematics            
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            else
                if nargin == 3
                    obj.Kinematics(q, dq, zeros([obj.n, 1], 'like', q));
                else
                    if nargin == 4
                        obj.Kinematics(q, dq, ddq);
                    end
                end
            end
            % Compute the term
            int_r_X_ddr_ = zeros(3, 1, 'like', obj.rGauss(1:3, 1));
            for i = 1:obj.NGaussPointsInt
                int_r_X_ddr_ = int_r_X_ddr_ + cross(obj.rGauss(:, i), obj.ddrGauss(:, i))*obj.LinearMassDensity(obj.GaussPoints(i))*obj.GaussWeights(i);
            end
        end
        
        % Integral of \cross(\dor{r}, \jacobian{r}{q})
        function int_dr_X_pv_r_ = int_dr_X_pv_r(obj, q, dq)
            if nargin == 2 % Update the kinematics
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            else
                if nargin == 3
                    obj.Kinematics(q, dq, zeros([obj.n, 1], 'like', q));
                end
            end
            % Compute the term
            int_dr_X_pv_r_ = zeros(obj.n, 3, 'like', obj.rGauss(1:3, 1));
            for i = 1:obj.NGaussPointsInt
                int_dr_X_pv_r_ = int_dr_X_pv_r_ + (skew(obj.drGauss(:, i))*obj.JdrGauss(1:3, 1:obj.n, i))'*obj.LinearMassDensity(obj.GaussPoints(i))*obj.GaussWeights(i);
            end
        end

        %Integral of \dot(\jacobian{r}{q}, \ddot{r})
        function int_pv_r_O_dd_r_ = int_pv_r_O_dd_r(obj, q, dq, ddq)
            if nargin == 2 % Update the kinematics
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            else
                if nargin == 3
                    obj.Kinematics(q, dq, zeros([obj.n, 1], 'like', q));
                else
                    if nargin == 4
                        obj.Kinematics(q, dq, ddq);
                    end
                end
            end
            % Compute the term
            int_pv_r_O_dd_r_ = zeros(obj.n, 1, 'like', obj.rGauss(1:3, 1));
            for i = 1:obj.NGaussPointsInt
                int_pv_r_O_dd_r_ = int_pv_r_O_dd_r_ + (obj.JdrGauss(1:3, 1:obj.n, i)'*obj.ddrGauss(:, i))*obj.LinearMassDensity(obj.GaussPoints(i))*obj.GaussWeights(i);
            end
        end

        % Integral of \dot(\dot{r}, \dot{r})
        % TODO: This is always zero, remove
        function int_dr_O_dr_ = int_dr_O_dr(obj, q, dq)
            int_dr_O_dr_ = zeros(1, 1, 'like', q);
        end
        
        % Jacobian of the integral of \dot{r}
        % TODO: This is always zero, remove
        function grad_int_dr_ = grad_int_dr(obj, q)
            grad_int_dr_ = zeros(obj.n, 3, 'like', q);
        end
        
        % Jacobian of the integral of \cross{r, \dot{r}}
        function grad_int_r_X_dr_ = grad_int_r_X_dr(obj, q)
            if nargin == 2% Update the DK
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end
            % Compute the term
            grad_int_r_X_dr_ = zeros(obj.n, 3, 'like', obj.rGauss(1:3, 1));
            for i = 1:obj.NGaussPointsInt
                grad_int_r_X_dr_ = grad_int_r_X_dr_ +  (skew(obj.rGauss(:, i))*obj.JdrGauss(1:3, 1:obj.n, i))'*obj.LinearMassDensity(obj.GaussPoints(i))*obj.GaussWeights(i);
            end
        end

        %Jacobian of the time derivative of the inertia
        function grad_J_ = grad_J(obj, q)
            if nargin == 2 % Update the DK
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end
            % Compute the term
            grad_J_ = zeros(3, 3, obj.n, 'like', obj.rGauss(1:3, 1));
            % Iterate over all DoFs
            for j = 1:obj.n
                % Iterate over all the Gaussian points
                for i = 1:obj.NGaussPointsInt
                    % i-th Gaussian point
                    ri                   = obj.rGauss(1:3, i);
                    % Gradient of i-th Gaussian point with respect to q(j)
                    dri_qj               = obj.JdrGauss(1:3, j, i);
                    grad_J_(1:3, 1:3, j) = grad_J_(1:3, 1:3, j) + (skew(dri_qj)'*skew(ri) + skew(ri)'*skew(dri_qj))*obj.GaussWeights(i)*obj.LinearMassDensity(obj.GaussPoints(i));
                end
            end
        end
        
        % Jacobian of the center of mass velocity
        function grad_v_com_ = grad_v_com(obj, q)
            if nargin == 2 %Update the DK
                obj.Kinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end
            grad_v_com_ = obj.grad_v_com_;
        end

        % Strain function
        function xi_ = xi(obj, q, s)
            xi_ = obj.StrainBasis(s)*q + obj.ReferenceStrain;
        end

        % Strain Jacobian with respect to q
        function Jxi_ = Jxi(obj, q, s)
            Jxi_ = obj.StrainBasis(s);
        end

        % First time derivative of the strain
        function dxi_ = dxi(obj, q, dq, s)
            dxi_ = obj.Jxi(q, s)*dq;
        end

        % Second time derivative of the strain
        function ddxi_ = ddxi(obj, q, dq, ddq, s)
            ddxi_ = obj.Jxi(q, s)*ddq;
        end

        % Generalized elastic force
        function K_ = K(obj, q)
            G   = obj.ShearModulus;
            E   = obj.YoungModulus;
            K_  = zeros(obj.n, 1, 'like', q);
            for i = 1:obj.NGaussPointsInt
                % Linear and angular stiffness
                K_l         = pi*obj.Radius(obj.GaussPoints(i))^2*diag([G, G, E]);
                BodyInertia = pi*obj.Radius(obj.GaussPoints(i))^4*[1/4, 1/4, 1/2];
                K_a         = diag([E*BodyInertia(1), E*BodyInertia(2), G*BodyInertia(3)]);
                % Body stiffness
                K_b         = blkdiag(K_a, K_l);
                % Update the elastic force
                xiGauss     = obj.xi(q, obj.GaussPoints(i));
                JxiGauss    = obj.Jxi(q, obj.GaussPoints(i));
                K_          = K_ + (JxiGauss'*K_b*(xiGauss - obj.ReferenceStrain))*obj.GaussWeights(i);
            end
        end
        % Generalized damping force
        function D_ = D(obj, q, dq)
            G   = obj.ShearModulus;
            E   = obj.YoungModulus;
            D_ = zeros(obj.n, 1, 'like', q);
            for i = 1:obj.NGaussPointsInt
                % Linear and angular damping
                D_l         = pi*obj.Radius(obj.GaussPoints(i))^2*diag([G, G, E]);
                BodyInertia = pi*obj.Radius(obj.GaussPoints(i))^4*[1/4, 1/4, 1/2];
                D_a         = diag([E*BodyInertia(1), E*BodyInertia(2), G*BodyInertia(3)]);
                % Body damping
                D_b         = blkdiag(D_a, D_l)*obj.MaterialDamping;
                % Update the damping force
                dxiGauss    = obj.dxi(q, dq, obj.GaussPoints(i));
                JxiGauss    = obj.Jxi(q, obj.GaussPoints(i));
                D_          = D_ + (JxiGauss'*D_b*dxiGauss)*obj.GaussWeights(i);
            end
        end
    end

    %% Private methods
    methods (Access = private)
        % Compute the kinematics
        function [g, omega, v, domega, dv, J_omega, J_v] = Kinematics(obj, q, dq, ddq)
            % Compute the direct and differential kinamtics at the tip
            [g, omega, v, domega, dv, J_omega, J_v] = obj.Kinematics_s(q, dq, ddq, obj.RestLength);
            % Update the position of Gauss points of the backbone in the tip frame
            RT          = g(1:3, 1:3)';
            dRT         = -RT*skew(obj.EtaGauss(1:3, end));
            ddRT        = -RT*skew(obj.dEtaGauss(1:3, end)) + ...
                          -dRT*skew(obj.EtaGauss(1:3, end));
            d           = g(1:3, 4);
            pGauss      = zeros(3, obj.NGaussPointsInt, 'like', q);
            dpGauss     = zeros(3, obj.NGaussPointsInt, 'like', q);
            JdpGauss    = zeros(3, obj.n, obj.NGaussPointsInt, 'like', q);
            ddpGauss    = zeros(3, obj.NGaussPointsInt, 'like', q);
            
            if coder.target("MATLAB")
                obj.rGauss              = zeros(3, obj.NGaussPointsInt, 'like', q);
                obj.drGauss             = zeros(3, obj.NGaussPointsInt, 'like', q);
                obj.JdrGauss            = zeros(3, obj.n, obj.NGaussPointsInt, 'like', q);
                obj.ddrGauss            = zeros(3, obj.NGaussPointsInt, 'like', q);
            end
            for i = 1:obj.NGaussPointsInt
                pGauss(:, i)    = RT*(obj.gGauss(1:3, 4, i) - d);
                dpGauss(:, i)   = (dRT*(obj.gGauss(1:3, 4, i) - obj.gGauss(1:3, 4, end)) + ...
                                         RT*(obj.EtaGauss(4:6, i)  - obj.EtaGauss(4:6, end)));
                JdpGauss(:, 1:obj.n, i) = RT*(obj.JEtaGauss(4:6, 1:obj.n, i) + skew(obj.gGauss(1:3, 4, i))*obj.JEtaGauss(1:3, 1:obj.n, end));
                ddpGauss(:, i)  = (ddRT*(obj.gGauss(1:3, 4, i) - obj.gGauss(1:3, 4, end)) + ...
                                         2*dRT*(obj.EtaGauss(4:6, i)  - obj.EtaGauss(4:6, end)) + ...
                                          RT*(obj.dEtaGauss(4:6, i)  - obj.dEtaGauss(4:6, end)));
            end
            %% Update the CoM position and its time derivatives
            % Compute the CoM position
            obj.p_com_  = zeros(3, 1, 'like', q);
            m       = obj.m();
            for i = 1:obj.NGaussPointsInt
                obj.p_com_  = obj.p_com_ + pGauss(:, i)*obj.GaussWeights(i)*obj.LinearMassDensity(obj.GaussPoints(i));
            end
            obj.p_com_ = (1/m)*obj.p_com_;
            % Compute the CoM velocity
            obj.v_com_rel_   = zeros(3, 1, 'like', q);
            for i = 1:obj.NGaussPointsInt
                obj.v_com_rel_  = obj.v_com_rel_ + obj.GaussWeights(i)*obj.LinearMassDensity(obj.GaussPoints(i))*dpGauss(:, i);
            end
            obj.v_com_rel_      = (1/m)*obj.v_com_rel_;
            % Compute the CoM acceleration
            obj.a_com_rel_ = zeros(3, 1, 'like', q);
            for i = 1:obj.NGaussPointsInt
                obj.a_com_rel_  = obj.a_com_rel_ + obj.GaussWeights(i)*obj.LinearMassDensity(obj.GaussPoints(i))*ddpGauss(:, i);
            end
            obj.a_com_rel_ = (1/m)*obj.a_com_rel_;
            %% Update the relative position vector
            J_p_comi = zeros(3, obj.n, 'like', q);
            for i = 1:obj.NGaussPointsInt
                J_p_comi = J_p_comi + obj.GaussWeights(i)*obj.LinearMassDensity(obj.GaussPoints(i))*JdpGauss(1:3, 1:obj.n, i);
            end
            J_p_comi        = (1/m)*J_p_comi;
            for i = 1:obj.NGaussPointsInt
                obj.rGauss(1:3, i)              = pGauss(1:3, i) - obj.p_com_;
                obj.drGauss(1:3, i)             = dpGauss(1:3, i) - obj.v_com_rel_;
                obj.JdrGauss(1:3, 1:obj.n, i)   = JdpGauss(1:3, 1:obj.n, i) - J_p_comi;
                obj.ddrGauss(1:3, i)            = ddpGauss(1:3, i) - obj.a_com_rel_;
            end
            %% Update the gradient of the CoM velocity in the local frame
            obj.grad_v_com_ = skew(RT'*obj.p_com_)*obj.JEtaGauss(1:3, 1:obj.n, end) - obj.JEtaGauss(4:6, 1:obj.n, end);
            % Gradient of the CoM in the base frame
            grad_v_com_i_1  = zeros(3, obj.n, 'like', q);
            for i = 1:obj.NGaussPointsInt
                grad_v_com_i_1 = grad_v_com_i_1 + obj.JEtaGauss(4:6, 1:obj.n, i)*obj.GaussWeights(i)*obj.LinearMassDensity(obj.GaussPoints(i));
            end
            grad_v_com_i_1  = (1/m)*grad_v_com_i_1;
            obj.grad_v_com_ = (RT*(obj.grad_v_com_ + grad_v_com_i_1))';
        end
        %% Compute all the kinematic quantities at point s
        function [g, omega, v, domega, dv, J_omega, J_v] = Kinematics_s(obj, q, dq, ddq, s)
            % Transformation matrix from s to base 
            g        = eye(4, 'like', q);
            % Angular velocity
            omega    = zeros(3, 1, 'like', q);
            % Linear velocity
            v        = zeros(3, 1, 'like', q);
            % Angular acceleration
            domega    = zeros(3, 1, 'like', q);
            % Linear velocity
            dv        = zeros(3, 1, 'like', q);
            % Jacobians
            J_omega   = zeros(3, obj.n, 'like', q);
            J_v       = zeros(3, obj.n, 'like', q);
            % Body velocity 
            Eta      = zeros(6, 1, 'like', q);
            % Jacobian of body velocity
            JEta     = zeros(6, obj.n, 'like', q);
            % Body acceleration
            dEta     = zeros(6, 1, 'like', q);
            if coder.target("MATLAB")
                obj.gGauss              = repmat(eye(4, 'like', q), [1, 1, obj.NGaussPointsInt]);
                obj.EtaGauss            = zeros(6, obj.NGaussPointsInt, 'like', q);
                obj.JEtaGauss           = zeros(6, obj.n, obj.NGaussPointsInt, 'like', q);
                obj.dEtaGauss           = zeros(6, obj.NGaussPointsInt, 'like', q);
            end
            % Check trivial position, i.e., base
            if s == 0
                return;
            end
            % Iteration variables initialization
            exitFlag= false;
            hTotal  = 0;
            S       = 0;
            % Loop over all the Gaussian points
            for i = 1:obj.NGaussPointsInt
                % Compute the distance h between two Gaussian points
                if i == 1
                    h   = obj.GaussPoints(i);
                else
                    h   = obj.GaussPoints(i)-obj.GaussPoints(i-1);
                end
                % Check if s is in the current computation window
                if s >= hTotal && s <= hTotal + h
                    if i == 1
                        h       = s;
                    else
                        h       = s-obj.GaussPoints(i-1);
                    end
                    exitFlag= true;
                end
                % Evaluate the strain, its Jacobian and time derivatives at the Zanna quadrature points
                xi_1    = obj.xi(q, S + h/2 - sqrt(3)*h/6);
                xi_2    = obj.xi(q, S + h/2 + sqrt(3)*h/6);
                Jxi_1   = obj.Jxi(q, S + h/2 - sqrt(3)*h/6);
                Jxi_2   = obj.Jxi(q, S + h/2 + sqrt(3)*h/6);
                dxi_1   = Jxi_1*dq;
                dxi_2   = Jxi_2*dq;
                ddxi_1  = obj.ddxi(q, dq, ddq, S + h/2 - sqrt(3)*h/6);
                ddxi_2  = obj.ddxi(q, dq, ddq, S + h/2 + sqrt(3)*h/6);
                % Compute the Magnus expansion of the strain
                Omega   = (h/2)*(xi_1 + xi_2) + (sqrt(3)*h^2)/12*ad(xi_1)*xi_2;
                OmegaHat= [skew(Omega(1:3)), Omega(4:6); zeros(1, 4)];
                % Compute first order time derivative of Omega and its Jacobian with respect to q
                PhiOmega= (h/2)*(Jxi_1 + Jxi_2) +...
                          (sqrt(3)*h^2)/12*(ad(xi_1)*Jxi_2 - ad(xi_2)*Jxi_1);
                dOmega  = PhiOmega*dq;
                % Compute second order time derivative of Omega
                ddOmega = (h/2)*(ddxi_1 + ddxi_2) + ....
                            (sqrt(3)*h^2)/12*(ad(ddxi_1)*xi_2 + 2*ad(dxi_1)*dxi_2 + ad(xi_1)*ddxi_2);
                % Compute the tanget map
                TOmega  = Tanexpmat(Omega);
                % Update the transformation matrix and derivatives
                g       = g*expmat(OmegaHat);
                dEta    = invAd(expmat(OmegaHat))*(dEta -ad(TOmega*dOmega)*Eta + TOmega*ddOmega);
                JEta    = invAd(expmat(OmegaHat))*(JEta + TOmega*PhiOmega);
                Eta     = JEta*dq;
                % Store the values for future use but expressed in the
                % global frame
                obj.gGauss(:, :, i)                 = g;
                R                                   = g(1:3, 1:3);
                omega                               = R*Eta(1:3);
                dR                                  = skew(omega)*R;
                obj.EtaGauss(1:6, i)                = blkdiag(R, R)*Eta;
                obj.JEtaGauss(1:6, 1:obj.n, i)      = blkdiag(R, R)*JEta;
                obj.dEtaGauss(1:6, i)               = blkdiag(dR, dR)*Eta + blkdiag(R, R)*dEta;
                if exitFlag
                    break;
                end
                % Prepare for the next iteration
                hTotal = hTotal + h;
                S      = obj.GaussPoints(i);
            end
            % NOTE: All the results must be rotated in the intertial frame
            % since the above quantities are computed in body coordinates
            % Compute the angular velocity and its Jacobian
            R       = g(1:3, 1:3);
            J_omega = R*JEta(1:3, 1:obj.n);
            omega   = R*Eta(1:3);
            dR      = skew(omega)*R;
            % Compute the lienar tip velocity and its Jacobian
            J_v     = R*JEta(4:6, 1:obj.n);
            v       = R*Eta(4:6);
            % Compute the accelerations
            domega  = dR*Eta(1:3) + R*dEta(1:3);
            dv      = dR*Eta(4:6) + R*dEta(4:6);
        end
    end
end

