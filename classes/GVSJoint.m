classdef GVSJoint < JointNew
    %GVSBODY Class representing a slender body modeled under the geometric
    %variable strain approach.

    methods (Abstract)
        %Strain basis as a function of the arc length
        StrainBasis(obj, s);
    end
    
    properties
        %Parameters are in this order:
        %L_0 = rest length
        %R   = radius
        %rho = mass density
        %E   = Young's modulus
        %Poi = Poisson's ratio
        %Eta = Material damping
        %GaussPoints = Number of Gaussian Points
        Parameters;
        RestLength  = 0;
        %Reference strain
        ReferenceStrain = [0;0;0;0;0;1];
    end

    %Private properties extracted from the parameters for ease of use
    properties (Access = private)
        Radius          = 0;
        MassDensity     = 0;
        YoungModulus    = 0;
        ShearModulus    = 0;
        PoissonRatio    = 0.5;
        MaterialDamping = 0;
        NGaussPoints    = 0;
        GaussPoints     = 0;
        GaussWeights    = 0;
        NGaussPointsInt = 0;
        %Properties used for the internal computation, all the quantities
        %are expressed in the global frame
        gGauss;
        EtaGauss;
        JEtaGauss;
        dEtaGauss;
        rGauss;
        drGauss;
        JdrGauss;
        ddrGauss;
        LinearMassDensity = 0;%Mass density aling the curvilinear abscissa
    end
    
    methods
        %% Class constructor
        function obj = GVSJoint(n, Parameters)
            if isrow(Parameters)
                Parameters = Parameters';
            end
            %Call the superclass constructor
            obj                     = obj@JointNew(n);
            %Set the parameters
            obj.Parameters          = Parameters;
            obj.RestLength          = Parameters(1);
            obj.Radius              = Parameters(2);
            obj.MassDensity         = Parameters(3);
            obj.LinearMassDensity   = pi*obj.Radius^2*obj.MassDensity;
            obj.YoungModulus        = Parameters(4);
            obj.PoissonRatio        = Parameters(5);
            obj.ShearModulus        = obj.YoungModulus/(2*(1+obj.PoissonRatio));
            obj.MaterialDamping     = Parameters(6);
            obj.NGaussPoints        = Parameters(7);
            [obj.GaussPoints, obj.GaussWeights] = lgwt(obj.NGaussPoints, 0, obj.RestLength);
            obj.GaussPoints         = [obj.GaussPoints; obj.RestLength];
            obj.GaussWeights        = [obj.GaussWeights; 0];
            obj.NGaussPointsInt     = obj.NGaussPoints + 1;
            %obj.n                   = size(obj.StrainBasis(0), 2);
            % if isstring(StrainBasis) || ischar(StrainBasis)
            %     obj.StrainBasis         = str2func(StrainBasis);
            % else
            %     obj.StrainBasis         = StrainBasis;
            % end
            obj.gGauss              = repmat(eye(4)     , [1, 1, obj.NGaussPointsInt]);
            obj.EtaGauss            = zeros(6, obj.NGaussPointsInt);
            obj.JEtaGauss           = zeros(6, obj.n, obj.NGaussPointsInt);
            obj.dEtaGauss           = zeros(6, obj.NGaussPointsInt);
            obj.rGauss              = zeros(3, obj.NGaussPointsInt);
            obj.drGauss             = zeros(3, obj.NGaussPointsInt);
            obj.JdrGauss            = zeros(3, obj.n, obj.NGaussPointsInt);
            obj.ddrGauss            = zeros(3, obj.NGaussPointsInt);
        end
        %% Methods implementation
        %Overload the update body method
        function updateJoint(obj, q, dq, ddq)
            %Update the kinematics
            [T_, omega_, v_, domega_, dv_, J_omega, J_v] = obj.Kinematics(q, dq, ddq);
            obj.T_                  = T_;
            obj.v_rel_              = v_;
            obj.omega_rel_          = omega_;
            obj.a_rel_              = dv_;
            obj.domega_rel_         = domega_;
            RT                      = T_(1:3, 1:3)';
            obj.v_par_              = RT*J_v;
            obj.omega_par_          = RT*J_omega;
            % disp("[GVS Joint] Update");
            % obj.T_
            % obj.v_rel_
            % obj.omega_rel_
            % obj.a_rel_
            % obj.domega_rel_
            % obj.v_par_
            % obj.omega_par_
            % disp("[GVS Joint] End update");
        end
        %% Transformation from base to body tip
        function T_ = T(obj, q)
            T_ = obj.T_s(q, obj.RestLength);
        end
        %% Transformation matrix from base to s
        function Ts_ = T_s(obj, q, s)
            [Ts_, ~, ~, ~, ~, ~, ~] = obj.Kinematics_s(q, zeros(obj.n, 1, 'like', q), zeros(obj.n, 1, 'like', q), s);
        end
        %% Relative velocity of the body tip
        function v_rel_ = v_rel(obj, q, dq)
            [~, ~, v_rel_, ~, ~, ~, ~] = obj.Kinematics(q, dq, zeros(obj.n, 1, 'like', q));
        end
        %% Relative angular velocity
        function omega_rel_ = omega_rel(obj, q, dq)
            [~, omega_rel_, ~, ~, ~, ~, ~] = obj.Kinematics(q, dq, zeros(obj.n, 1, 'like', q));
        end
        %% Relative linear acceleration of the center of mass
        function a_rel_ = a_rel(obj, q, dq, ddq)
            [~, ~, ~, ~, a_rel_, ~, ~] = obj.Kinematics(q, dq, ddq);
        end
        %% Relative angular acceleration
        function domega_rel_ = domega_rel(obj, q, dq, ddq)
            [~, ~, ~, domega_rel_, ~, ~, ~] = obj.Kinematics(q, dq, ddq);
        end
        %% Jacobian of the linear velocity of the tip with respect to q
        function v_par_ = v_par(obj, q)
            [g, ~, ~, ~, ~, ~, v_par_] = obj.Kinematics(q, zeros(obj.n, 1, 'like', q), zeros(obj.n, 1, 'like', q));
            v_par_ = g(1:3, 1:3)'*v_par_;
        end
        %% Jacobian of the angular velocity with respect to q
        function omega_par_ = omega_par(obj, q)
            [g, ~, ~, ~, ~, omega_par_, ~] = obj.Kinematics(q, zeros(obj.n, 1, 'like', q), zeros(obj.n, 1, 'like', q));
            omega_par_ = g(1:3, 1:3)'*omega_par_;
        end
        %% Strain function
        function xi_ = xi(obj, q, s)
            xi_ = obj.StrainBasis(s)*q + obj.ReferenceStrain;
        end
        %% Strain Jacobian with respect to q
        function Jxi_ = Jxi(obj, q, s)
            Jxi_ = obj.StrainBasis(s);
        end
        %% First time derivative of the strain
        function dxi_ = dxi(obj, q, dq, s)
            dxi_ = obj.Jxi(q, s)*dq;
        end
        %% Second time derivative of the strain
        function ddxi_ = ddxi(obj, q, dq, ddq, s)
            ddxi_ = obj.Jxi(q, s)*ddq;
        end
    end

    %% Private methods
    methods (Access = private)
        %% Compute the kinematics
        function [g, omega, v, domega, dv, J_omega, J_v] = Kinematics(obj, q, dq, ddq)
            %% Compute the direct and differential kinamtics at the tip
            [g, omega, v, domega, dv, J_omega, J_v] = obj.Kinematics_s(q, dq, ddq, obj.RestLength);
        end
        %% Compute all the kinematic quantities at point s
        function [g, omega, v, domega, dv, J_omega, J_v] = Kinematics_s(obj, q, dq, ddq, s)
            %Preallocate the outputs
            %Transformation matrix from s to base 
            g        = eye(4, 'like', q);
            %Angular velocity
            omega    = zeros(3, 1, 'like', q);
            %Linear velocity
            v        = zeros(3, 1, 'like', q);
            %Angular acceleration
            domega    = zeros(3, 1, 'like', q);
            %Linear velocity
            dv        = zeros(3, 1, 'like', q);
            %Jacobians
            J_omega   = zeros(3, obj.n, 'like', q);
            J_v       = zeros(3, obj.n, 'like', q);
            %Body velocity 
            Eta      = zeros(6, 1, 'like', q);
            %Jacobian of body velocity
            JEta     = zeros(6, obj.n, 'like', q);
            %Body acceleration
            dEta     = zeros(6, 1, 'like', q);
            %Check trivial position, i.e., base
            if s == 0
                return;
            end
            %Iteration variables initialization
            exitFlag= false;
            hTotal  = 0;
            S       = 0;
            %Loop over all the Gaussian points
            for i = 1:obj.NGaussPointsInt
                %Compute the distance h between two Gaussian points
                if i == 1
                    h   = obj.GaussPoints(i);
                else
                    h   = obj.GaussPoints(i)-obj.GaussPoints(i-1);
                end
                %Check if s is in the current computation window
                if s >= hTotal && s <= hTotal + h
                    if i == 1
                        h       = s;
                    else
                        h       = s-obj.GaussPoints(i-1);
                    end
                    exitFlag= true;
                end
                %Evaluate the strain, its Jacobian and time derivatives at the Zanna quadrature points
                xi_1    = obj.xi(q, S + h/2 - sqrt(3)*h/6);
                xi_2    = obj.xi(q, S + h/2 + sqrt(3)*h/6);
                Jxi_1   = obj.Jxi(q, S + h/2 - sqrt(3)*h/6);
                Jxi_2   = obj.Jxi(q, S + h/2 + sqrt(3)*h/6);
                dxi_1   = Jxi_1*dq;
                dxi_2   = Jxi_2*dq;
                ddxi_1  = obj.ddxi(q, dq, ddq, S + h/2 - sqrt(3)*h/6);
                ddxi_2  = obj.ddxi(q, dq, ddq, S + h/2 + sqrt(3)*h/6);
                %Compute the Magnus expansion of the strain
                Omega   = (h/2)*(xi_1 + xi_2) + (sqrt(3)*h^2)/12*ad(xi_1)*xi_2;
                OmegaHat= [skew(Omega(1:3)), Omega(4:6); zeros(1, 4)];
                %Compute first order time derivative of Omega and its Jacobian with respect to q
                PhiOmega= (h/2)*(Jxi_1 + Jxi_2) +...
                          (sqrt(3)*h^2)/12*(ad(xi_1)*Jxi_2 - ad(xi_2)*Jxi_1);
                dOmega  = PhiOmega*dq;
                %Compute second order time derivative of Omega
                ddOmega = (h/2)*(ddxi_1 + ddxi_2) + ....
                            (sqrt(3)*h^2)/12*(ad(ddxi_1)*xi_2 + 2*ad(dxi_1)*dxi_2 + ad(xi_1)*ddxi_2);
                %Compute the tanget map
                TOmega  = Tanexpmat(Omega);
                %Update the transformation matrix and derivatives
                g       = g*expmat(OmegaHat);
                dEta    = invAd(expmat(OmegaHat))*(dEta -ad(TOmega*dOmega)*Eta + TOmega*ddOmega);
                JEta    = invAd(expmat(OmegaHat))*(JEta + TOmega*PhiOmega);
                Eta     = JEta*dq;
                %Store the values for future use but expressed in the
                %global frame
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
                %Prepare for the next iteration
                hTotal = hTotal + h;
                S      = obj.GaussPoints(i);
            end
            %NOTE: All the results must be rotated in the intertial frame
            %since the above quantities are computed in body coordinates
            %Compute the angular velocity and its Jacobian
            R       = g(1:3, 1:3);
            J_omega = R*JEta(1:3, 1:obj.n);
            omega   = R*Eta(1:3);
            dR      = skew(omega)*R;
            %Compute the lienar tip velocity and its Jacobian
            J_v     = R*JEta(4:6, 1:obj.n);
            v       = R*Eta(4:6);
            %Compute the accelerations
            domega  = dR*Eta(1:3) + R*dEta(1:3);
            dv      = dR*Eta(4:6) + R*dEta(4:6);
        end
    end
end

