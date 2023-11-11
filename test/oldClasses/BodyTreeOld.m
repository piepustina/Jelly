classdef BodyTreeOld < handle
    %BODYTREEOLD Class representing an open tree of rigid and soft bodies
    %#codegen

    %!!!!!!!!!!!!!!!!!!!!!! THIS CLASS HAS TO BE REMOVED

    properties (Constant)
        %Define the maximum number of bodies and joints (required for code generation)
        MaxBodiesNumber = 20;
    end

    properties (Access = public)
        %Joints of the tree
        Joints;
        %Bodies of the tree
        Bodies;
        %Total number of bodies in the tree
        N_B = 0;
        %Total number of DoFs
        n    = 0;
        %Gravity field in the base frame
        g    = [0; 0; -9.81];
        %Orientation of the base with respect to the world frame
        T0   = eye(4);
        %Mass condition numer
        MassConditionNumber = 0;
    end

    properties (Access = public)
        %Threshold for the mass matrix and Coriolis term
        Q_THSD_M     = 0;
        %Reference value for q for the Taylor approximation of M and C
        Q_THSD_M_REF = 0;
        %Threshold for the gravity vector
        Q_THSD_G     = 0;
        %Reference value for q for the Taylor approximation of G
        Q_THSD_G_REF = 0;
    end

    methods
        function obj = BodyTreeOld(Joints, Bodies, varargin)
            %BODYTREE Construct a BodyTree consisting of at most N joints and
            %bodies whose type is specified by the string arrays joints and
            %bodies
            %Transform the input to handle the code generation issues
            %
            computeConfigurationLimits = false;
            switch nargin
                case 2
                    computeConfigurationLimits = true;
                case 6
                    obj.Q_THSD_M = varargin{1};
                    obj.Q_THSD_M_REF = varargin{2};
                    obj.Q_THSD_G = varargin{3};
                    obj.Q_THSD_G_REF = varargin{4};
                otherwise
                    error("Incorrect number of arguments.");
            end
            obj.N_B = 0;
            obj.n   = 0;
            l_B     = length(Bodies);
            for i = 1:BodyTree.MaxBodiesNumber
                if i <= l_B
                    if ~isnumeric(Bodies{i})
                        obj.n = obj.n + Bodies{i}.n;
                        obj.N_B = obj.N_B + 1;
                    end
                end
            end
            obj.Joints = cell(BodyTree.MaxBodiesNumber, 1);
            obj.Bodies = cell(BodyTree.MaxBodiesNumber, 1);
            obj.Joints = Joints;
            obj.Bodies = Bodies;
            %Compute the configuration tresholds if not provided as
            %arguments
            if computeConfigurationLimits
                obj.Q_THSD_M = zeros(obj.n, 1);
                obj.Q_THSD_M_REF = zeros(obj.n, 1);
                obj.Q_THSD_G = zeros(obj.n, 1);
                obj.Q_THSD_G_REF = zeros(obj.n, 1);
                obj.ComputeConfigurationThresholds();
            end
        end
        
        %Evaluate the joints and bodies in the current configuration
        function obj = TreeUpdate(obj, q, dq, ddq)
            %Initialize variables
            k_j = 1;
            k_b = 1;
            l_B = length(obj.Bodies);
            for i = 1:BodyTree.MaxBodiesNumber
                if i <= l_B
                    if isnumeric(obj.Bodies{i})
                        continue;
                    end
                    %JOINT UPDATE
                    %Get the DOFs associated with the joint
                    n_j   = obj.Joints{i}.n;
                    %Compute the last index of q associated with the current
                    %joint
                    k_j_1 = k_j + n_j - 1;
                    %Compute (q, dq, ddq) associated with the joint
                    q_j   = q(k_j:k_j_1);
                    dq_j  = dq(k_j:k_j_1);
                    ddq_j = ddq(k_j:k_j_1);
                    %Update the joint data
                    obj.Joints{i}.updateJoint(q_j, dq_j, ddq_j);
                    %Update the joint index for the next iteration
                    k_j = k_j_1 + 1;
                    
                    %BODY UPDATE
                    %Get the DOFs associated with the body
                    n_b   = obj.Bodies{i}.n;
                    %Compute the last index of q associated with the current
                    %body
                    k_b_1 = k_b + n_b - 1;
                    %Compute (q, dq, ddq) associated with the body
                    q_b   = q(k_b:k_b_1);
                    dq_b  = dq(k_b:k_b_1);
                    ddq_b = ddq(k_b:k_b_1);
                    %Update the body data
                    obj.Bodies{i}.updateBody(q_b, dq_b, ddq_b);
                    %Update the body index for the next iteration
                    k_b = k_b_1 + 1;
                end
            end
        end
        
        %Implementation of the Inverse Dynamics using Kane's equations
        %type = 'double' if q is numeric, 'sym' if q is symbolic
        %NOTE: The inverse dynamics does not account for external and
        %generalized forces, except for gravity.
        function tau = InverseDynamics(obj, q, dq, ddq, type)
            %Update the BodyTree
            obj = obj.TreeUpdate(q, dq, ddq);
            %Run the ID
            tau = obj.Kane_aux(obj.g, type);
        end

        function tau = Kane_aux(obj, g, type)
            %KANE_AUX Given a obj model compute the generalized force that realizes ddq in
            %the state (q,dq)
            
            %Define auxiliary variables
            N_B_ = obj.N_B;
            %Define the output
            %tau = cell(N_B, 1);
            tau = zeros(obj.n, 1, type);
            %3 x N_B matrix whose i-th column stores the linear velocity of Body i ( origin of frame {S_i} )
            v      = zeros(3, N_B_, type);
            %3 x N_B matrix whose i-th column stores the angular velocity of Body i
            omega  = zeros(3, N_B_, type);
            %3 x N_B matrix whose i-th column stores the linear acceleration of Body i ( origin of frame {S_i} )
            a      = zeros(3, N_B_, type);
            %3 x N_B matrix whose i-th column stores the angular acceleration of Body i
            domega = zeros(3, N_B_, type);
            %3 x N_B matrix whose i-th column stores the linear acceleration of CoM_i
            a_com  = zeros(3, N_B_, type);
            %3 x N_B matrix whose i-th column stores vector Gamma_i
            Gamma  = zeros(3, N_B_, type);
            %3 x N_B matrix whose i-th column stores vector Omega_i
            Omega  = zeros(3, N_B_, type);
            %3 x N_B matrix whose i-th column stores vector M_i
            M      = zeros(3, N_B_, type);
            %3 x N_B matrix whose i-th column stores vector N_i
            N      = zeros(3, N_B_, type);
            %Auxiliary vector to represent the velocity of the preceding vector
            v_i_1       = zeros(3, 1, type);
            omega_i_1   = zeros(3, 1, type);
            a_i_1       = -g;
            domega_i_1  = zeros(3, 1, type);
            M_i_1  = zeros(3, 1, type);
            N_i_1  = zeros(3, 1, type);
            %********************************************************
            %********************* Forward step *********************
            %********************************************************
            %for i = 1:N_B_ %Iterative over all the bodies
            l_B = length(obj.Bodies);
            for i = 1:BodyTree.MaxBodiesNumber %Iterative over all the bodies
                if i <= l_B
                    if isnumeric(obj.Bodies{i}) || isnumeric(obj.Joints{i})
                        continue;
                    end
                    %Step 1
                    R_i_T        = real(obj.Joints{i}.T_(1:3, 1:3)');
                    t_i          = real(obj.Joints{i}.T_(1:3, 4));
                    v_rel_i      = real(obj.Joints{i}.v_rel_);
                    omega_rel_i  = real(obj.Joints{i}.omega_rel_);
                    dv_rel_i     = real(obj.Joints{i}.a_rel_);
                    domega_rel_i = real(obj.Joints{i}.domega_rel_);
                
                    p_com_i      = real(obj.Bodies{i}.p_com_);
                    v_com_rel_i  = real(obj.Bodies{i}.v_com_rel_);
                    a_com_rel_i  = real(obj.Bodies{i}.a_com_rel_);
                    
                    v(:, i)      = real(R_i_T*(v_i_1 + cross(omega_i_1, t_i) + v_rel_i));
                    omega(:, i)  = real(R_i_T*(omega_i_1 + omega_rel_i));
                    a(:, i)      = real(R_i_T*(a_i_1 + ...
                                          cross(domega_i_1, t_i) + ...
                                          cross(omega_i_1, cross(omega_i_1, t_i) + v_rel_i) + ...
                                          cross(omega_i_1, v_rel_i) + ...
                                          dv_rel_i));
                    domega(:, i) = real(R_i_T*(domega_i_1 + cross(omega_i_1, omega_rel_i) + domega_rel_i));
                    a_com(:, i) = real(a(:, i) + cross(domega(:, i), p_com_i) + cross(omega(:, i), cross(omega(:, i), p_com_i) + v_com_rel_i) +...
                                    cross(omega(:, i), v_com_rel_i) + a_com_rel_i);
                    
                    %Step 2
                    Gamma(:, i) =real(a_com(:, i)*obj.Bodies{i}.m_ +...
                                    2*cross(omega(:, i), obj.Bodies{i}.int_dr_) +...
                                    obj.Bodies{i}.int_ddr_);
                    Omega(:, i) = real(obj.Bodies{i}.I_*domega(:, i) + cross(omega(:, i), obj.Bodies{i}.I_*omega(:, i)) +...
                                    obj.Bodies{i}.J_*omega(:, i)+...
                                    cross(omega(:, i), obj.Bodies{i}.int_r_X_dr_) + ...
                                    obj.Bodies{i}.int_r_X_ddr_);
                    
                    %Update the iteration variables
                    v_i_1 = v(:, i);
                    omega_i_1 = omega(:, i);
                    a_i_1 = a(:, i);
                    domega_i_1 = domega(:, i);
                end
            end
            
            %********************************************************
            %********************* Backward step ********************
            %********************************************************
            %Index for tau vector
            idx_tau = obj.n;
            %for i = N_B_:-1:1
            for i = BodyTree.MaxBodiesNumber:-1:1
                if i <= l_B
                    if isnumeric(obj.Bodies{i}) || isnumeric(obj.Joints{i})
                        continue;
                    end
                    %Step 3
                    if i == N_B_
                        M(:, i) = real(Gamma(:, i));
                        N(:, i) = real(Omega(:, i) + cross(obj.Bodies{i}.p_com_, Gamma(:, i)));
                    else
                        if ~isnumeric(obj.Joints{i+1})
                            R_i_1   = real(obj.Joints{i+1}.T_(1:3, 1:3));
                            t_i_1   = real(obj.Joints{i+1}.T_(1:3, 4));
                            M(:, i) = real(Gamma(:, i) + R_i_1*M_i_1);
                            N(:, i) = real(Omega(:, i) + cross(obj.Bodies{i}.p_com_, Gamma(:, i)) + R_i_1*N_i_1 + cross(t_i_1, R_i_1*M_i_1));
                        end
                    end
                    
                    %Step 4
                    tau(idx_tau - obj.Bodies{i}.n + 1:idx_tau) = real(obj.Bodies{i}.grad_int_dr_*a_com(:, i) +...
                                                                     obj.Bodies{i}.grad_int_r_X_dr_*domega(:, i) +...
                                                                     arrayfun(@(j) -(1/2)*omega(:, i)'*obj.Bodies{i}.grad_J_(:, :, j)*omega(:, i), (1:size(obj.Bodies{i}.grad_J_, 3))')+...
                                                                     2*obj.Bodies{i}.int_dr_X_pv_r_*omega(:, i) +...
                                                                     obj.Bodies{i}.int_pv_r_O_dd_r_ +...
                                                                     obj.Bodies{i}.grad_v_com_*Gamma(:, i) +...
                                                                     (obj.Joints{i}.v_par_')*M(:, i) +...
                                                                     (obj.Joints{i}.omega_par_')*N(:, i));                   
                    %Update the iteration variables
                    M_i_1 = M(:, i);
                    N_i_1 = N(:, i);
                    idx_tau = idx_tau - obj.Bodies{i}.n;
                end
            end
        end
        
        %Forward dynamics in state space form
        function dx = StateSpaceForwardDynamics(obj, t, x, u, type)
            switch nargin
                case 1 || 2
                    x = zeros(2*obj.n, 1);
                    u = zeros(obj.n, 1);
                    type = 'double';
                case 3
                    u = zeros(obj.n, 1);
                    type = 'double';
                case 4
                    type = 'dobule';
            end
            q  = x(1:obj.n);
            dq = x(obj.n+1:end);
            dx = [dq; obj.ForwardDynamics(q, dq, u, type)];
        end
        
        %Implements the Forward Dynamics using Kane's equations. Cost is
        %cubic in the worst case.
        function ddq = ForwardDynamics(obj, q, dq, u, type)
            %Compute the mass matrix
            M = obj.MassMatrix(q, type);
            %Improve the condition number of M
            c = obj.MassConditionNumber;
            M = M + c*eye(size(M));
            %Coriolis and gravitational terms
            %Less efficient but handles better equilibria because different
            %values are used
            CG = obj.ApparentForce(q, dq, type) + obj.GravityForce(q, type);
            %Overall forces
            f  =  -CG - obj.K(q, type) - obj.D(q, dq, type) + u;
            %Perform scaling to improve simulation accuracy
            %Compute the acceleration
            ddq = pinv(M)*f;
            %ddq = M\f;
        end

        %Compute the generalized elastic force in the current
        %configuration.
        %To implement a custom value, just overload this method.
        function Kq = K(obj, q, type)
            %Update the tree only if q is passed as argument
            switch nargin
                case 1
                    type = 'double';
                case 2
                    obj.TreeUpdate(q, zeros(obj.n, 1), zeros(obj.n, 1));
                    type = 'double';
                case 3
                    obj.TreeUpdate(q, zeros(obj.n, 1), zeros(obj.n, 1));
            end
            Kq = zeros(obj.n, 1, type);
            k_b = obj.n;
            l_B = length(obj.Bodies);
            for i = BodyTree.MaxBodiesNumber:-1:1 %Iterative over all the bodies
                if i <= l_B
                    if isnumeric(obj.Bodies{i}) || isnumeric(obj.Joints{i})
                        continue;
                    end
                    k_b_1 = k_b - obj.Bodies{i}.n + 1;
                    %Kq(k_b_1:k_b) = real(obj.Bodies{i}.K_);
                    Kq(k_b_1:k_b) = obj.Bodies{i}.K_;
                    %Prepare for the next iteration
                    k_b = k_b_1 - 1;
                end
            end
        end

        %Compute the generalized damping force in the current
        %configuration.
        %To implement a custom value, just overload this method.
        function Dq = D(obj, q, dq, type)
            %Update the tree only if q and dq are passed as arguments
            switch nargin
                case 1
                    type = 'double';
                case 2
                    type = 'double';
                case 3
                    obj.TreeUpdate(q, dq, zeros(obj.n, 1));
                    type = 'double';
                case 4
                    obj.TreeUpdate(q, dq, zeros(obj.n, 1));
            end
            Dq = zeros(obj.n, 1, type);
            k_b = obj.n;
            l_B = length(obj.Bodies);
            for i = BodyTree.MaxBodiesNumber:-1:1 %Iterative over all the bodies discarding the fake ones
                if i <= l_B
                    if isnumeric(obj.Bodies{i}) || isnumeric(obj.Joints{i})
                        continue;
                    end
                    k_b_1 = k_b - obj.Bodies{i}.n + 1;
                    %Dq(k_b_1:k_b) = real(obj.Bodies{i}.D_);
                    Dq(k_b_1:k_b) = obj.Bodies{i}.D_;
                    %Prepare for the next iteration
                    k_b = k_b_1 - 1;
                end
            end
        end

        %Mass matrix.
        %NOTE: Current implementation has cost O(n^2)!
        function M = MassMatrix(obj, q, type)
            %Store value of gravity and set gravity to zero to compute the
            %mass matrix
            g_ = obj.g;
            obj.g = zeros(3, 1);
            M = zeros(obj.n, obj.n, type);
            Zeron = zeros(obj.n, 1);
            Id   = eye(obj.n);
            for i = 1:obj.n
                %M(:, i) = obj.RobustInverseDynamics(q, Zeron, Id(:, i), type, obj.Q_THSD_M, obj.Q_THSD_M_REF);
                M(:, i) = obj.InverseDynamics(q, Zeron, Id(:, i), type);
            end
            %Symmetrize the result because of possible numerical
            %approximations
            M = 1/2*(M + M');
            %Restrore gravity and compute the generalized forces
            obj.g = g_;
        end

        %Gravity force
        function G = GravityForce(obj, q, type)
            Zeron = zeros(obj.n, 1);
            G     = obj.RobustInverseDynamics(q, Zeron, Zeron, type, obj.Q_THSD_G, obj.Q_THSD_G_REF);
        end

        function C = ApparentForce(obj, q, dq, type)
            %Set gravity force to zero
            g_ = obj.g;
            obj.g = zeros(3, 1);
            %Compute the apparent forces
            C = obj.RobustInverseDynamics(q, dq, zeros(obj.n, 1), type, obj.Q_THSD_M, obj.Q_THSD_M_REF);
            %Restore back gravity
            obj.g = g_;
        end
        
        
        %Computes the apparent (Coriolis/Centrifugal) matrix. We use the
        %approach of [Kawasaki et al., IEEE T-RA 1996]. The cost is
        %quadratic.
        function C = ApparentMatrix(obj, q, dq, type)
            %Set gravity force to zero
            g_ = obj.g;
            obj.g = zeros(3, 1);
            %Preallocate the apparent matrix
            C = zeros(obj.n, obj.n, type);
            In= eye(obj.n);
            %Compute the apparent forces
            for i = 1:obj.n
                C(:, i) = 1/2*( obj.RobustInverseDynamics(q, dq + In(:, i), zeros(obj.n, 1), type, obj.Q_THSD_M, obj.Q_THSD_M_REF) ...
                               -obj.RobustInverseDynamics(q, dq           , zeros(obj.n, 1), type, obj.Q_THSD_M, obj.Q_THSD_M_REF) ...
                               -obj.RobustInverseDynamics(q, In(:, i)     , zeros(obj.n, 1), type, obj.Q_THSD_M, obj.Q_THSD_M_REF));
            end
            %Restore back gravity
            obj.g = g_;
        end
        
        %Compute the equilibrium configuration from q0 and given the
        %actuation tau
        function [q_eq, f_val] = EquilibriumConfiguration(obj, q0, tau)
            %options = optimoptions('fmincon', 'FunctionTolerance', 1e-4);
            %q_eq = fsolve(@(q) obj.EquilibriumEquation(q, tau), q0, options);
            [q_eq, f_val] = fsolve(@(q) obj.EquilibriumEquation(q, tau), q0);
        end
    end

    methods (Access = protected)
        %Check if we need the computation of a limit
        function [q_approx, approximated] = ApproxQ(~, q, thsd, q_thsd)
            q_approx = q;
            approximated = false;
            for i = 1:length(q)
                if abs(q(i)) <= thsd(i)
                    approximated = true;
                    if q(i) ~= 0
                        q_approx(i) = q_thsd(i)*sign(q(i));
                    else
                        q_approx(i) = q_thsd(i);
                    end
                end
            end
        end

        %A robustified version of the inverse dynamics that uses
        %a Taylor expansion of order 1 when the configuration is close to a
        %limit value
        function tau = RobustInverseDynamics(obj, q, dq, ddq, type, thsd, q_thsd)
            
            if strcmp(type, 'sym')
                %Call directly the inverse dynamics
                tau = obj.InverseDynamics(q, dq, ddq, type);
            else
                [Q, approximated] = obj.ApproxQ(q, thsd, q_thsd);
                %If no approximation of q is required just use the
                %procedure, otherwise Taylor approximate
                if ~approximated
                    tau = obj.InverseDynamics(q, dq, ddq, type);
                else
                    %Compute tau at Q
                    tauQ              = obj.InverseDynamics(Q, dq, ddq, type);
                    %Compute the Jacobian of the inverse dynamics in Q
                    epsilon           = min(q_thsd);%Do a step that is safe for the computation of the inverse dynamics
                    JtauQ             = obj.InverseDynamicsJacobian(Q, dq, ddq, epsilon, tauQ);
                    tau               = tauQ + JtauQ*(q - Q);
                end
            end
        end
        
        %Implements the equlibrium equation
        function eq = EquilibriumEquation(obj, q, tau)
            eq = obj.GravityForce(q, 'double') + obj.K(q, 'double') - tau;
        end
        

        %Numerical computation of the jacobian of the inverse dynamics
        %w.r.t. q only. The Jacobian is evaluated in q and epsilon is the
        %step size for the computation of the Jacobian. When epsilon is 0,
        %the function returns a matrix of zeros by default. 
        %tauq is the value of the inverse dynamics in q so that we do not
        %have to repear the computation.
        function jac = InverseDynamicsJacobian(obj, q, dq, ddq, epsilon, tauq)
        %Output initialization
        jac = zeros(obj.n, obj.n);
        if epsilon == 0
            return;
        end
        
        %Inverse of the step size
        epsilon_inv = 1/epsilon;
        
        % Do perturbation with a forward approach
        for i = 1 : obj.n
            %Perturb the configuration
            q_plus     = q;
            q_plus(i)  =  q(i) + sign(q(i))*epsilon;
            %Run the ID in q_plus
            tau_plus = obj.InverseDynamics(q_plus, dq, ddq, 'double');
            %Approximate the derivative of the inverse dynamics w.r.t. q(i)
            jac(:, i) = epsilon_inv * sign(q(i))*(tau_plus - tauq);
        end
        end
    
        %Compute the thresholds required for the computations of the limits
        %around the origin
        function ComputeConfigurationThresholds(obj)
            %Hyperparameters
            q_ref = 0.2;%reference value to use for the mass matrix and gravity vector
            %Compute the limits for M and C starting from 0
            %The evaluation is made with steps of 10^{-1} starting from
            %zero up to 10^{-16}
            q          = q_ref*ones(obj.n, 1);
            %Compute the order of magnitude of the norm of M in the reference value q
            M_norm_ref = norm(obj.MassMatrix(q, 'double'));
            n_M_norm_ref = floor(log10(M_norm_ref));
            for i = 1:obj.n
                %Compute the reference value for q(i)
                q(i) = obj.Q_THSD_M_REF(i);
                for j = -16:-1
                    %As a measure of how good the approximation is we use
                    %the norm of M
                    M_q      = obj.MassMatrix(q, 'double');
                    M_norm   = norm(M_q);
                    lambda   = eig(M_q);
                    n_M_Norm = floor(log10(M_norm));
                    if n_M_Norm <= n_M_norm_ref && all(lambda > 0)
                        break;
                    else
                        q(i) = 10^(j);
                    end
                end
                obj.Q_THSD_M_REF(i) = q(i);
                %Compute the interval in which to use the approximation
                obj.Q_THSD_M(i)     = q(i);
                n_M_norm_tshd = n_M_Norm;
                for k=j:-1
                    M_norm   = norm(obj.MassMatrix(q, 'double'));
                    n_M_Norm = floor(log10(M_norm));
                    if n_M_Norm <= n_M_norm_tshd
                        break;
                    else
                        q(i) = 10^(k);
                    end
                end
                obj.Q_THSD_M(i)     = q(i);
                %Reset q
                q(i) = q_ref;
            end
            %Repeat the process for the gravitational force
            %Compute the order of magnitude of the norm of M in the reference value q
            G_norm_ref = norm(obj.GravityForce(q, 'double'));
            n_G_norm_ref = floor(log10(G_norm_ref));
            for i = 1:obj.n
                %Compute the reference value for q(i)
                q(i) = obj.Q_THSD_G_REF(i);
                for j = -16:-1
                    %As a measure of how good the approximation is we use
                    %the norm of M
                    G_norm   = norm(obj.GravityForce(q, 'double'));
                    n_G_Norm = floor(log10(G_norm));
                    if n_G_Norm <= n_G_norm_ref
                        break;
                    else
                        q(i) = 10^(j);
                    end
                end
                obj.Q_THSD_G_REF(i) = q(i);
                %Compute the interval in which to use the approximation
                obj.Q_THSD_G(i)     = q(i);
                n_G_norm_thsd = n_G_Norm;
                for k=j:-1
                    G_norm   = norm(obj.GravityForce(q, 'double'));
                    n_G_Norm = floor(log10(G_norm));
                    if n_G_Norm <= n_G_norm_thsd
                        break;
                    else
                        q(i) = 10^(k);
                    end
                end
                obj.Q_THSD_G(i)     = q(i);
                %Reset q
                q(i) = q_ref;
            end
        end
    end
end

