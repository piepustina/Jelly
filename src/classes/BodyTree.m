classdef BodyTree < handle
    %Class modeling a mechanical system of generic joints and bodies serially interconnected. 

    %#codegen

    % TODO: Convert Joints and Bodies from cell arrays to simple arrays,
    % apparently this works for code generation as well. 

    properties (Constant)
        % Maximum number of bodies and joints, required for code generation.
        MaxBodiesNumber = 20;
    end

    properties (Access = public)
        %Joints of the kinematic tree organized in a cell array.
        Joints;
        %Bodies of the kinematic tree organized in a cell array.
        Bodies;
        %Total number of joints/bodies in the kinematic tree.
        N_B = 0;
        %Total number of DoF of the kinematic tree.
        n    = 0;
        %Gravity force in the base frame.
        g    = [0; 0; -9.81];
        %Orientation of the base frame with respect to the world frame.
        T0   = eye(4);
        
        MassConditionNumber = 0;%TODO: Remove
    end
    
    properties (Access = private)
        % Joints can be treated as bodies with no inertial parameters. Thus,
        % the augments the state by considering also the joints as bodies.
        BodiesInternal;
        N_B_Internal = 0;
    end
    
    methods
        function obj = BodyTree(Joints, Bodies)
            %Construct a BodyTree consisting of joints and
            %bodies.
            %
            %Args:
                %    Joints ({:class:`Joint`}): Cell array of joints of the tree
                %    Bodies ({:class:`Body`}) : Cell array of bodies of the tree
            
            % Compute the number of bodies
            obj.N_B = 0;
            obj.n   = 0;
            l_B     = length(Bodies);
            for i = 1:BodyTree.MaxBodiesNumber
                if i <= l_B
                    if ~isnumeric(Bodies{i})
                        obj.n = obj.n + Bodies{i}.n + Joints{i}.n;
                        obj.N_B = obj.N_B + 1;
                    end
                end
            end
            % Assign the joints and bodies
            obj.Joints = cell(BodyTree.MaxBodiesNumber, 1);
            obj.Bodies = cell(BodyTree.MaxBodiesNumber, 1);
            obj.Joints = Joints;
            obj.Bodies = Bodies;
        
            % Allocate the augmeneted bodies for code generation
            obj.N_B_Internal    = 2*obj.N_B;
            obj.BodiesInternal  = cell(2*BodyTree.MaxBodiesNumber, 1);
            j                   = 1;
            if coder.target("MATLAB")
                for i = 1:2:obj.N_B_Internal
                    obj.BodiesInternal{i}   = obj.Joints{j};
                    obj.BodiesInternal{i+1} = obj.Bodies{j};
                    j = j + 1;
                end
            else
                for i = 1:2:2*BodyTree.MaxBodiesNumber
                   obj.BodiesInternal{i}   = obj.Joints{j};
                   obj.BodiesInternal{i+1} = obj.Bodies{j};
                   j = j + 1;
                end
            end
        end
        
        function obj = TreeUpdate(obj, q, dq, ddq)
            %Update the state of the BodyTree. 
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of configuration variables

            % Check that the vectors q, dq and ddq are columns.
            if ~iscolumn(q)
                q = q';
            end
            if ~iscolumn(dq)
                dq = dq';
            end
            if ~iscolumn(ddq)
                ddq = ddq';
            end
            
            % Initialize variables
            k_i = 1;
            l_B = length(obj.Bodies);
            for i = 1:BodyTree.MaxBodiesNumber
                if i <= l_B
                    if isnumeric(obj.Bodies{i})
                        continue;
                    end
                    % JOINT UPDATE
                    % Get the DOFs associated with the joint
                    n_j   = obj.Joints{i}.n;
                    if n_j ~= 0
                        % Compute the last index of q associated with the current
                        % joint
                        k_i_1 = k_i + n_j - 1;
                        % Compute (q, dq, ddq) associated with the joint
                        q_j   = q(k_i:k_i_1);
                        dq_j  = dq(k_i:k_i_1);
                        ddq_j = ddq(k_i:k_i_1);
                        % Update the joint data
                        obj.Joints{i}.Update(q_j, dq_j, ddq_j);
                        % Update the joint index for the next iteration
                        k_i = k_i_1 + 1;
                    elseobj.N_B;
                        obj.Joints{i}.Update([], [], []);
                    end
                    
                    % BODY UPDATE
                    % Get the DOFs associated with the body
                    n_b   = obj.Bodies{i}.n;
                    if n_b ~= 0
                        % Compute the last index of q associated with the current
                        % body
                        k_i_1 = k_i + n_b - 1;
                        % Compute (q, dq, ddq) associated with the body
                        q_b   = q(k_i:k_i_1);
                        dq_b  = dq(k_i:k_i_1);
                        ddq_b = ddq(k_i:k_i_1);
                        % Update the body data
                        obj.Bodies{i}.Update(q_b, dq_b, ddq_b);
                        % Update the body index for the next iteration
                        k_i = k_i_1 + 1;
                    else
                        obj.Bodies{i}.Update([], [], []);
                    end
                end
            end
        end
        
        function tau = InverseDynamics(obj, q, dq, ddq)
            % Evaluate the inverse dynamics using the Generalized Inverse Dynamics (GID) algorithm.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of configuration variables

            % Check that the input vectors are in column format.
            if ~iscolumn(q)
                q = q';
            end
            if ~iscolumn(dq)
                dq = dq';
            end
            if ~iscolumn(ddq)
                ddq = ddq';
            end
            % Update the BodyTree
            obj = obj.TreeUpdate(q, dq, ddq);
            % Run the GID algorithm,  q is only passed to
            % retreive its type.
            tau = obj.Kane_aux(obj.g, q);
        end
        
        function dx = StateSpaceForwardDynamics(obj, t, x, u)
            %Evaluate the forward dynamics in state space form.
            %
            %Args:
            %   t (double): Time
            %   x ([double], [sym]): State variables consisting of configuration variables and time derivative of the configuration variables
            %   u ([double], [sym]): Generalized actuation force

            switch nargin
                case 1 || 2
                    x = zeros(2*obj.n, 1);
                    u = zeros(obj.n, 1);
                case 3
                    u = zeros(obj.n, 1, "like", x);
            end
            q  = x(1:obj.n);
            dq = x(obj.n+1:end);
            dx = [dq; obj.ForwardDynamics(q, dq, u)];
        end
        
        function ddq = ForwardDynamics(obj, q, dq, tau)
            %Forward dynamics through the inverse dynamics algorithm.
            %
            %Args:
            %   q   ([double], [sym]): Configuration variables
            %   dq  ([double], [sym]): First-order time derivative of the configuration variables
            %   tau ([double], [sym]): Generalized actuation force
            
            % Compute the mass matrix
            M = obj.MassMatrix(q);
            % Improve the condition number of M
            c = obj.MassConditionNumber;
            M = M + c*eye(size(M), "like", q);
            % Coriolis and gravitational terms
            % Less efficient but handles better equilibria because different
            % values are used
            CG = obj.ApparentForce(q, dq) + obj.GravityForce(q);
            % Overall forces
            f  =  -CG - obj.K(q) - obj.D(q, dq) + tau;
            % Perform scaling to improve simulation accuracy
            % Compute the acceleration
            ddq = M\f;
        end

        function Kq = K(obj, q)
            %Evaluate the generalized elastic force.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables

            % Update the tree only if q is passed as argument
            switch nargin
                case 1
                    q = zeros(obj.n, 1);
                case 2
                    obj.TreeUpdate(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q));
            end
            Kq = zeros(obj.n, 1, "like", q);
            k_b = obj.n;
            l_B = obj.N_B_Internal;
            for i = 2*BodyTree.MaxBodiesNumber:-1:1 % Iterative over all the augmeneted bodies
                if i <= l_B
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    if obj.BodiesInternal{i}.n ~= 0
                        k_b_1 = k_b - obj.BodiesInternal{i}.n + 1;
                        Kq(k_b_1:k_b) = obj.BodiesInternal{i}.K_;
                        % Prepare for the next iteration
                        k_b = k_b_1 - 1;
                    end
                end
            end
        end

        function Dq = D(obj, q, dq)
            %Evaluate the generalized damping force.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables

            % Update the tree only if q and dq are passed as arguments
            switch nargin
                case 1
                    q = zeros(obj.n, 1);
                case 2
                    obj.TreeUpdate(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q));
                case 3
                    obj.TreeUpdate(q, dq, zeros(obj.n, 1, "like", q));
            end
            Dq = zeros(obj.n, 1, "like", q);
            k_b = obj.n;
            l_B = obj.N_B_Internal;
            for i = 2*BodyTree.MaxBodiesNumber:-1:1 % Iterative over all the bodies discarding the fake ones
                if i <= l_B
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    if obj.BodiesInternal{i}.n ~= 0
                        k_b_1 = k_b - obj.BodiesInternal{i}.n + 1;
                        Dq(k_b_1:k_b) = obj.BodiesInternal{i}.D_;
                        % Prepare for the next iteration
                        k_b = k_b_1 - 1;
                    end
                end
            end
        end

        function M = MassMatrix(obj, q)
            %Evaluate the mass matrix.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
            
            % Store value of gravity and set gravity to zero to compute the
            % mass matrix
            g_      = obj.g;
            obj.g   = zeros(3, 1, "like", g_);
            M       = zeros(obj.n, obj.n, "like", q);
            Zeron   = zeros(obj.n, 1, "like", q);
            Id      = eye(obj.n, "like", q);
            for i = 1:obj.n
                M(:, i) = obj.InverseDynamics(q, Zeron, Id(:, i));
            end
            % Symmetrize the result because of possible numerical
            % approximations
            M       = 1/2*(M + M');
            % Restrore gravity and compute the generalized forces
            obj.g   = g_;
        end

        function G = GravityForce(obj, q)
            %Evaluate the generalized gravitational force.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables

            Zeron = zeros(obj.n, 1, "like", q);
            G     = obj.InverseDynamics(q, Zeron, Zeron);
        end

        function C = ApparentForce(obj, q, dq)
            %Evaluate the generalized apparent force.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables

            % Set gravity force to zero
            g_ = obj.g;
            obj.g = zeros(3, 1, "like", g_);
            % Compute the apparent forces
            C = obj.InverseDynamics(q, dq, zeros(obj.n, 1, "like", q));
            % Restore back gravity
            obj.g = g_;
        end
        
        
        function C = ApparentMatrix(obj, q, dq)
            %Evaluate the apparent matrix.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables

            % We use the
            % approach of [Kawasaki et al., IEEE T-RA 1996]. The cost is
            % quadratic.
            % 
            % Set gravity force to zero
            g_      = obj.g;
            obj.g   = zeros(3, 1, "like", g_);
            % Preallocate the apparent matrix
            C       = zeros(obj.n, obj.n, "like", q);
            In      = eye(obj.n, "like", q);
            Zeron   = zeros(obj.n, 1, "like", q);
            % Compute the apparent forces
            for i = 1:obj.n
                C(:, i) = 1/2*( obj.InverseDynamics(q, dq + In(:, i), Zeron) ...
                               -obj.InverseDynamics(q, dq           , Zeron) ...
                               -obj.InverseDynamics(q, In(:, i)     , Zeron));
            end
            % Restore back gravity
            obj.g = g_;
        end
        
        function [q_eq, f_val] = EquilibriumConfiguration(obj, q0, tau)
            %Find an equilibrium configuration by solving numerically the equilibrium equations.
            %
            %Args:
                %    q0  ([double]): Initial guess for the equlibrium
                %    tau ([double]): Generalized actuation force

            % options = optimoptions('fmincon', 'FunctionTolerance', 1e-4);
            % q_eq = fsolve(@(q) obj.EquilibriumEquation(q, tau), q0, options);
            [q_eq, f_val] = fsolve(@(q) obj.EquilibriumEquation(q, tau), q0);
        end
        

        function E = Energy(obj, q, dq, q_ref)
            %Evaluate the system energy. It is assumed that the elastic
            %force, if any is linear.
            %
            %Args:
                %    q      ([double], [sym]): Configuration variables
                %    dq     ([double], [sym]): First-order time derivative of configuration variables
                %    q_ref  ([double], [sym]): Reference configuration variables for the elastic and gravitational energy
            
            % q_ref is the reference configuration for the computation of
            % the elastic and gravitational energy. If it is not provided as
            % argument, 0 is assumed.
            switch nargin
                case 2
                    dq   = zeros(obj.n, 1, "like", q);
                    q_ref = zeros(obj.n, 1, "like", q);
                case 3
                    q_ref = zeros(obj.n, 1, "like", q);
            end

            % Kinetic energy
            E_kinetic = 1/2*dq'*obj.MassMatrix(q)*dq;
            % Potential energy
            E_elastic = 1/2*(q - q_ref)'*obj.K(q);
            if ~isnumeric(q)
                E_gravity = potential(obj.GravityForce(q), q, q_ref);
            else
                E_gravity = integral(@(s) (q - q_ref)'*obj.GravityForce(q_ref + s.*(q - q_ref)), 0, 1, 'ArrayValued', true);
            end
            % Overall energy
            E = E_kinetic + E_elastic + E_gravity;
        end

        function T = DirectKinematics(obj, q, idx)
            %Evaluate the direct kinematics.
            %
            %Args:
            %   q   ([double], [sym]): Configuration variables
            %   idx           ([int]): Ordered array of body indexes indicating the bodies for which the jacobian has to be evaluated
            %Return:
            %   ([double], [sym]): Homogeneous transformation matrices for the body specified by idx.   

            switch nargin
                case 2
                    idx = linspace(1, obj.N_B, obj.N_B);
            end

            % Modify the index to account for the fact that internally the joints are modeled as bodies
            idx = 2*idx;

            % Check that the input vectors are in column format.
            if ~iscolumn(q)
                q = q';
            end

            % T is matrix of vertically stacked 4x4 transformation matrices
            T   = repmat(eye(4, 'like', q), length(idx), 1);
            
            % Update the state of the kinematic tree
            Zeron   = zeros(obj.n, 1, "like", q);
            obj     = obj.TreeUpdate(q, Zeron, Zeron);
            
            % Variables initialization
            T_i         = obj.T0;
            j           = 1;
            lastIdx     = idx(end);
            N_B_        = obj.N_B_Internal;
            % Compute the direct kinematics
            for i = 1:2*BodyTree.MaxBodiesNumber % Iterative over all the augmented bodies
                if i <= N_B_
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end

                    T_i  = T_i*obj.BodiesInternal{i}.T_;
                    
                    if i == idx(j)
                        T(1+4*(j-1):4*j, 1:4)   = T_i;
                        j                       = j + 1;
                    end
    
                    % If we have reached the last body, break the loop
                    if i == lastIdx
                        break;
                    end
                end
            end
        end

        function [q, converged, e] = InverseKinematics(obj, T, idx, q0, N, task_flags)
            %Evaluate the inverse kinematics numerically using a Newton-Rapson iteration scheme.
            %
            %Args:
            %   T   ([double, double])      : Target transformation matrices vertically stacked
            %   q0  ([double, double])      : Initial guess, the default value is q0 = zeros(n, 1)
            %   idx ([double])              : Vector of indexes identifing the bodies in the chain for which T has to be found
            %   N   (double)                : Maximum number of iterations
            %   task_flags ([double, bool]) : Vector of flags specifying for each indexed body what components of the task vector should be considered
            %Return:
            %   {[double], [sym]}: Homogeneous transformation matrices for each body.   
            
            % Default values
            DefaultN            = 4;   %Number of Newton iterations
            AngularErrorThsd    = 1e-3;%Threshold in the Newton scheme for the angular velocity
            LinearErrorThsd     = 1e-3;%Threshold in the Newton scheme for the linear velocity
            
            switch nargin
                case 2
                    idx         = linspace(1, obj.N_B, obj.N_B);
                    q0          = zeros(obj.n, 1);
                    N           = DefaultN;
                    task_flags  = ones(obj.N_B*6, 1);
                case 3
                    q0  = zeros(obj.n, 1);
                    N   = DefaultN;
                    task_flags  = ones(obj.N_B*6, 1);
                case 4
                    N   = DefaultN;
                    task_flags  = ones(obj.N_B*6, 1);
                case 5
                    task_flags  = ones(obj.N_B*6, 1);
            end

            % Store useful variables
            idxLength = length(idx);

            % Convert the task flags to an index vector
            task_flags_l = logical(task_flags == 0);
            
            % Check that the dimensions of T and idx are consistent. 
            if floor(size(T, 1)/4) ~= idxLength
                error("The size of T and idx is not consistent.");
            end
            % Check that the dimensions of idx and task_flags are consistent
            if 6*idxLength ~= length(task_flags)
                error("The length of the body indexes is not consitent with the length of the task flags.");
            end

            % Preallocate the output for code generation
            q         = q0;
            converged = 0;

            % Allocate iteration variables
            e         = Inf*ones(6*idxLength, 1);
            e_thsd    = ones(idxLength, 1);%If contains all zeros, the configuration satisfies all the constraints

            % Run the Newton algorithm as given in Linch and Park, Modern Robotics
            for i = 1:N
                % Evaluate the direct kinematics in the current configuration
                T_q         = obj.DirectKinematics(q, idx);
                % Evaluate the body Jacobian in the current configuration
                J_q         = obj.BodyJacobian(q, idx);
                
                % Iterate over all the requried bodies
                for j = 1:idxLength
                    % Evaluate the desired configuration in the current frame
                    T_qd_j = invTransformation(T_q(1+4*(j-1):4*j, 1:4))*T(1+4*(j-1):4*j, 1:4);
                    
                    % Compute the error
                    e_j                                 = skew4_inv(logmat(T_qd_j));
                    % Set to zero the components for which the error does not matter
                    e_j(task_flags_l(1+6*(j-1):6*j))    = 0;
                    % Store the error to use in the update phase
                    e(1+6*(j-1):6*j)                    = e_j;
                    % Check if the error is below the threshold
                    if (norm(e_j(1:3)) <= AngularErrorThsd) && (norm(e_j(4:6)) <= LinearErrorThsd)
                        e_thsd(j) = 0;
                    end
                end

                % Check if the error is small enough on all the channels and in case exit, iteration converged
                if ~any(e_thsd)
                    converged = 1;
                    break;
                else
                    % Update the configuration
                    q = q + J_q\e;
                end
            end

            % Display a warning if the convergence was not achieved
            if ~converged
                warning("Convergence not achived. The error norm is: " + sprintf("%.4f", norm(e)));
            end

        end

        function J = BodyJacobian(obj, q, idx)
            %Evaluate the body Jacobian of the i-th body. If i is not specified, the method returns the body Jacobian of each body of the chain.
            %
            %Args:
            %   q   ([double], [sym]): Configuration variables
            %   idx           ([int]): Array of body indexes indicating the bodies for which the jacobian has to be evaluated
            %Return:
            %   {[double], [sym]}: length(idx)*6 x n body Jacobian with angular and linear components for each body specified by idx

            switch nargin
                case 1
                    q   = zeros(obj.n, 1, 'like', q);
                    idx = linspace(1, obj.N_B, obj.N_B);
                case 2
                    idx = linspace(1, obj.N_B, obj.N_B);
            end

            % Modify the index to account for the fact that internally the joints are modeled as bodies
            idx = 2*idx;

            % Check that the input vectors are in column format.
            if ~iscolumn(q)
                q = q';
            end

            % Update the BodyTree
            Zeron   = zeros(obj.n, 1, 'like', q);
            obj     = obj.TreeUpdate(q, Zeron, Zeron);

            % Define auxiliary variables
            % The joints are treated as massless bodies thus we augment the
            % body dimensions
            N_B_ = obj.N_B_Internal;
            
            % Define the output
            J      = zeros(length(idx)*6, obj.n, 'like', q);
            
            % Uncomment below for testing purposes
            % v      = zeros(3, 1, "like", q);
            % omega  = zeros(3, 1, "like", q);
            
            % Auxiliary variables for the iteration
            J_i        = zeros(6, obj.n, 'like', q);
            q_idx      = 1;
            j          = 1;
            lastIdx    = idx(end);
            
            % Compute the body Jacobian recursively
            for i = 1:2*BodyTree.MaxBodiesNumber % Iterative over all the augmented bodies
                if i <= N_B_
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    % Retreive the required information from the body
                    R_i_T           = real(obj.BodiesInternal{i}.T_(1:3, 1:3)');
                    t_i             = real(obj.BodiesInternal{i}.T_(1:3, 4));
                    v_par_i         = real(obj.BodiesInternal{i}.v_par_);
                    omega_par_i     = real(obj.BodiesInternal{i}.omega_par_);

                    % Update the variables
                    nBody                       = obj.BodiesInternal{i}.n;
                    if q_idx ~= 1
                        J_i(4:6, 1:q_idx-1)       = real(R_i_T*(J_i(4:6, 1:q_idx-1) - skew(t_i)*J_i(1:3, 1:q_idx-1)));
                        J_i(1:3, 1:q_idx-1)       = real(R_i_T*J_i(1:3, 1:q_idx-1));
                    end
                    J_i(1:6, q_idx:q_idx+nBody-1) = [omega_par_i; v_par_i];
                    
                    % Uncomment below for testing purposes
                    %v_rel_i         = real(obj.BodiesInternal{i}.v_rel_);
                    %omega_rel_i     = real(obj.BodiesInternal{i}.omega_rel_);
                    %v               = real(R_i_T*(v + cross(omega, t_i) + v_rel_i));
                    %omega           = real(R_i_T*(omega + omega_rel_i));

                    % Save the value if i hits the current body value
                    if i == idx(j)
                        J(1 + 6*(j-1):6*j, 1:obj.n) = J_i;
                        % Update the index for the body
                        j               = j + 1;
                    end

                    % Stop the iteration if the last body has been reached
                    if i == lastIdx
                        break;
                    end

                    % Prepare for the next iteration
                    q_idx           = q_idx + nBody;
                    
                end
            end
            
            % Uncomment below for testing purposes
            % disp("v");v
            % disp("omega");omega
            
        end
    end

    methods (Access = protected)
        
        function eq = EquilibriumEquation(obj, q, tau)
            %Evaluate the equilibrium equation.
            %
            %Args:
                %    q      ([double], [sym]): Configuration variables
                %    tau    ([double], [sym]): Generalized actuation force
            eq = obj.GravityForce(q) + obj.K(q) - tau;
        end
        
    end

    methods (Access = private)
        function tau = Kane_aux(obj, g, q_type)
            %Compute the generalized inverse dynamics using the recursion
            %formulas.
            %Args:
                %    g      ([double], [sym]): Gravitational force vector expressed in the base frame
                %    q_type ([double], [sym]): Dummy variable used to retrieve the type of the configuration variables
            
            % Define auxiliary variables
            % The joints are treated as massless bodies thus we augment the
            % body dimensions
            N_B_ = obj.N_B_Internal;
            
            % Define the output
            tau = zeros(obj.n, 1, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the linear velocity of Body i ( origin of frame {S_i} )
            v      = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the angular velocity of Body i
            omega  = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the linear acceleration of Body i ( origin of frame {S_i} )
            a      = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the angular acceleration of Body i
            domega = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the linear acceleration of CoM_i
            a_com  = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores vector Gamma_i
            Gamma  = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores vector Omega_i
            Omega  = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores vector M_i
            M      = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores vector N_i
            N      = zeros(3, N_B_, "like", q_type);
            % Auxiliary vector to represent the velocity of the preceding vector
            v_i_1       = zeros(3, 1, "like", q_type);
            omega_i_1   = zeros(3, 1, "like", q_type);
            a_i_1       = -g;
            domega_i_1  = zeros(3, 1, "like", q_type);
            M_i_1  = zeros(3, 1, "like", q_type);
            N_i_1  = zeros(3, 1, "like", q_type);
            % ********************************************************
            % ********************* Forward step *********************
            % ********************************************************
            for i = 1:2*BodyTree.MaxBodiesNumber % Iterative over all the augmented bodies
                if i <= N_B_
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    % Step 1
                    R_i_T        = real(obj.BodiesInternal{i}.T_(1:3, 1:3)');
                    t_i          = real(obj.BodiesInternal{i}.T_(1:3, 4));
                    v_rel_i      = real(obj.BodiesInternal{i}.v_rel_);
                    omega_rel_i  = real(obj.BodiesInternal{i}.omega_rel_);
                    dv_rel_i     = real(obj.BodiesInternal{i}.a_rel_);
                    domega_rel_i = real(obj.BodiesInternal{i}.domega_rel_);
                
                    p_com_i      = real(obj.BodiesInternal{i}.p_com_);
                    v_com_rel_i  = real(obj.BodiesInternal{i}.v_com_rel_);
                    a_com_rel_i  = real(obj.BodiesInternal{i}.a_com_rel_);
                    
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
                    
                    % Step 2
                    Gamma(:, i) =real(a_com(:, i)*obj.BodiesInternal{i}.m_ +...
                                    2*cross(omega(:, i), obj.BodiesInternal{i}.int_dr_) +...
                                    obj.BodiesInternal{i}.int_ddr_);
                    Omega(:, i) = real(obj.BodiesInternal{i}.I_*domega(:, i) + cross(omega(:, i), obj.BodiesInternal{i}.I_*omega(:, i)) +...
                                    obj.BodiesInternal{i}.J_*omega(:, i)+...
                                    cross(omega(:, i), obj.BodiesInternal{i}.int_r_X_dr_) + ...
                                    obj.BodiesInternal{i}.int_r_X_ddr_);
                    
                    % Update the iteration variables
                    v_i_1 = v(:, i);
                    omega_i_1 = omega(:, i);
                    a_i_1 = a(:, i);
                    domega_i_1 = domega(:, i);
                end
            end
            
            % ********************************************************
            % ********************* Backward step ********************
            % ********************************************************
            % Index for tau vector
            idx_tau = obj.n;
            for i = 2*BodyTree.MaxBodiesNumber:-1:1
                if i <= N_B_
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    % Step 3
                    if i == N_B_
                        M(:, i) = real(Gamma(:, i));
                        N(:, i) = real(Omega(:, i) + cross(obj.BodiesInternal{i}.p_com_, Gamma(:, i)));
                    else
                        if ~isnumeric(obj.BodiesInternal{i+1})
                            R_i_1   = real(obj.BodiesInternal{i+1}.T_(1:3, 1:3));
                            t_i_1   = real(obj.BodiesInternal{i+1}.T_(1:3, 4));
                            M(:, i) = real(Gamma(:, i) + R_i_1*M_i_1);
                            N(:, i) = real(Omega(:, i) + cross(obj.BodiesInternal{i}.p_com_, Gamma(:, i)) + R_i_1*N_i_1 + cross(t_i_1, R_i_1*M_i_1));
                        end
                    end
                    
                    % Step 4
                    if obj.BodiesInternal{i}.n ~= 0
                    tau(idx_tau - obj.BodiesInternal{i}.n + 1:idx_tau) = real(obj.BodiesInternal{i}.grad_int_dr_*a_com(:, i) +...
                                                                     obj.BodiesInternal{i}.grad_int_r_X_dr_*domega(:, i) +...
                                                                     arrayfun(@(j) -(1/2)*omega(:, i)'*obj.BodiesInternal{i}.grad_J_(:, :, j)*omega(:, i), (1:size(obj.BodiesInternal{i}.grad_J_, 3))')+...
                                                                     2*obj.BodiesInternal{i}.int_dr_X_pv_r_*omega(:, i) +...
                                                                     obj.BodiesInternal{i}.int_pv_r_O_dd_r_ +...
                                                                     obj.BodiesInternal{i}.grad_v_com_*Gamma(:, i) +...
                                                                     (obj.BodiesInternal{i}.v_par_')*M(:, i) +...
                                                                     (obj.BodiesInternal{i}.omega_par_')*N(:, i));                   
                    end
                    % Update the iteration variables
                    M_i_1 = M(:, i);
                    N_i_1 = N(:, i);
                    idx_tau = idx_tau - obj.BodiesInternal{i}.n;
                end
            end
        end
    end
end

