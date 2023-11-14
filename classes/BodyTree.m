classdef BodyTree < handle
    %Class representing an open tree of joints and bodies serially connected. 
    
    %#codegen

    %TODO: Convert Joints and Bodies from cell arrays to simple arrays,
    %apparently this works for code generation as well. 

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
    
    properties (Access = private)
        %Joints can be treated as bodies with no inertial parameters. Thus,
        %the augments the state by considering also the joints as bodies.
        BodiesInternal;
        N_B_Internal = 0;
    end
    
    methods
        %Class constructor
        function obj = BodyTree(Joints, Bodies)
            %Construct a BodyTree consisting of joints and
            %bodies.
            %
            %Args:
                %    Joints (Joint): Cell array of joints
                %    Bodies (Body): Cell array of bodies
                %               
            
            %Compute the number of bodies
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
            %Assign the joints and bodies
            obj.Joints = cell(BodyTree.MaxBodiesNumber, 1);
            obj.Bodies = cell(BodyTree.MaxBodiesNumber, 1);
            obj.Joints = Joints;
            obj.Bodies = Bodies;
            %Allocate the augmeneted bodies
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
        
        %Evaluate the joints and bodies in the current configuration
        function obj = TreeUpdate(obj, q, dq, ddq)
            %Update the state of the tree. 
            
            %Initialize variables
            k_i = 1;
            l_B = length(obj.Bodies);
            for i = 1:BodyTree.MaxBodiesNumber
                if i <= l_B
                    if isnumeric(obj.Bodies{i})
                        continue;
                    end
                    %JOINT UPDATE
                    %Get the DOFs associated with the joint
                    n_j   = obj.Joints{i}.n;
                    if n_j ~= 0
                        %Compute the last index of q associated with the current
                        %joint
                        k_i_1 = k_i + n_j - 1;
                        %Compute (q, dq, ddq) associated with the joint
                        q_j   = q(k_i:k_i_1);
                        dq_j  = dq(k_i:k_i_1);
                        ddq_j = ddq(k_i:k_i_1);
                        %Update the joint data
                        obj.Joints{i}.Update(q_j, dq_j, ddq_j);
                        %Update the joint index for the next iteration
                        k_i = k_i_1 + 1;
                    else
                        obj.Joints{i}.Update([], [], []);
                    end
                    
                    %BODY UPDATE
                    %Get the DOFs associated with the body
                    n_b   = obj.Bodies{i}.n;
                    if n_b ~= 0
                        %Compute the last index of q associated with the current
                        %body
                        k_i_1 = k_i + n_b - 1;
                        %Compute (q, dq, ddq) associated with the body
                        q_b   = q(k_i:k_i_1);
                        dq_b  = dq(k_i:k_i_1);
                        ddq_b = ddq(k_i:k_i_1);
                        %Update the body data
                        obj.Bodies{i}.Update(q_b, dq_b, ddq_b);
                        %Update the body index for the next iteration
                        k_i = k_i_1 + 1;
                    else
                        obj.Bodies{i}.Update([], [], []);
                    end
                end
            end
        end
        
        %Implementation of the Inverse Dynamics using Kane's equations
        %NOTE: The inverse dynamics does not account for external and
        %generalized forces, except for gravity.
        function tau = InverseDynamics(obj, q, dq, ddq)
            %Update the BodyTree
            obj = obj.TreeUpdate(q, dq, ddq);
            %Run the ID with the Kane algorithm. q is only passed to
            %retreive its type.
            tau = obj.Kane_aux(obj.g, q);
        end

        function tau = Kane_aux(obj, g, q_type)
            %KANE_AUX Given a obj model compute the generalized force that realizes ddq in
            %the state (q,dq)
            
            %Define auxiliary variables
            %The joints are treated as massless bodies thus we augment the
            %body dimensions
            N_B_ = obj.N_B_Internal;
            
            %Define the output
            tau = zeros(obj.n, 1, "like", q_type);
            %3 x N_B matrix whose i-th column stores the linear velocity of Body i ( origin of frame {S_i} )
            v      = zeros(3, N_B_, "like", q_type);
            %3 x N_B matrix whose i-th column stores the angular velocity of Body i
            omega  = zeros(3, N_B_, "like", q_type);
            %3 x N_B matrix whose i-th column stores the linear acceleration of Body i ( origin of frame {S_i} )
            a      = zeros(3, N_B_, "like", q_type);
            %3 x N_B matrix whose i-th column stores the angular acceleration of Body i
            domega = zeros(3, N_B_, "like", q_type);
            %3 x N_B matrix whose i-th column stores the linear acceleration of CoM_i
            a_com  = zeros(3, N_B_, "like", q_type);
            %3 x N_B matrix whose i-th column stores vector Gamma_i
            Gamma  = zeros(3, N_B_, "like", q_type);
            %3 x N_B matrix whose i-th column stores vector Omega_i
            Omega  = zeros(3, N_B_, "like", q_type);
            %3 x N_B matrix whose i-th column stores vector M_i
            M      = zeros(3, N_B_, "like", q_type);
            %3 x N_B matrix whose i-th column stores vector N_i
            N      = zeros(3, N_B_, "like", q_type);
            %Auxiliary vector to represent the velocity of the preceding vector
            v_i_1       = zeros(3, 1, "like", q_type);
            omega_i_1   = zeros(3, 1, "like", q_type);
            a_i_1       = -g;
            domega_i_1  = zeros(3, 1, "like", q_type);
            M_i_1  = zeros(3, 1, "like", q_type);
            N_i_1  = zeros(3, 1, "like", q_type);
            %********************************************************
            %********************* Forward step *********************
            %********************************************************
            l_B = length(obj.BodiesInternal);
            for i = 1:2*BodyTree.MaxBodiesNumber %Iterative over all the augmented bodies
                if i <= l_B
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    %Step 1
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
                    
                    %Step 2
                    Gamma(:, i) =real(a_com(:, i)*obj.BodiesInternal{i}.m_ +...
                                    2*cross(omega(:, i), obj.BodiesInternal{i}.int_dr_) +...
                                    obj.BodiesInternal{i}.int_ddr_);
                    Omega(:, i) = real(obj.BodiesInternal{i}.I_*domega(:, i) + cross(omega(:, i), obj.BodiesInternal{i}.I_*omega(:, i)) +...
                                    obj.BodiesInternal{i}.J_*omega(:, i)+...
                                    cross(omega(:, i), obj.BodiesInternal{i}.int_r_X_dr_) + ...
                                    obj.BodiesInternal{i}.int_r_X_ddr_);
                    
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
            for i = 2*BodyTree.MaxBodiesNumber:-1:1
                if i <= l_B
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    %Step 3
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
                    
                    %Step 4
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
                    %Update the iteration variables
                    M_i_1 = M(:, i);
                    N_i_1 = N(:, i);
                    idx_tau = idx_tau - obj.BodiesInternal{i}.n;
                end
            end
        end
        
        %Forward dynamics in state space form
        function dx = StateSpaceForwardDynamics(obj, ~, x, u)
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
        
        %Implements the Forward Dynamics using Kane's equations. Cost is
        %cubic in the worst case.
        function ddq = ForwardDynamics(obj, q, dq, u)
            %Compute the mass matrix
            M = obj.MassMatrix(q);
            %Improve the condition number of M
            c = obj.MassConditionNumber;
            M = M + c*eye(size(M), "like", q);
            %Coriolis and gravitational terms
            %Less efficient but handles better equilibria because different
            %values are used
            CG = obj.ApparentForce(q, dq) + obj.GravityForce(q);
            %Overall forces
            f  =  -CG - obj.K(q) - obj.D(q, dq) + u;
            %Perform scaling to improve simulation accuracy
            %Compute the acceleration
            ddq = pinv(M)*f;
            %ddq = M\f;
        end

        %Compute the generalized elastic force in the current
        %configuration.
        %To implement a custom value, just overload this method.
        function Kq = K(obj, q)
            %Update the tree only if q is passed as argument
            switch nargin
                case 1
                    q = zeros(obj.n, 1);
                case 2
                    obj.TreeUpdate(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q));
            end
            Kq = zeros(obj.n, 1, "like", q);
            k_b = obj.n;
            l_B = obj.N_B_Internal;
            for i = 2*BodyTree.MaxBodiesNumber:-1:1 %Iterative over all the augmeneted bodies
                if i <= l_B
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    if obj.BodiesInternal{i}.n ~= 0
                        k_b_1 = k_b - obj.BodiesInternal{i}.n + 1;
                        Kq(k_b_1:k_b) = obj.BodiesInternal{i}.K_;
                        %Prepare for the next iteration
                        k_b = k_b_1 - 1;
                    end
                end
            end
        end

        %Compute the generalized damping force in the current
        %configuration.
        %To implement a custom value, just overload this method.
        function Dq = D(obj, q, dq)
            %Update the tree only if q and dq are passed as arguments
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
            for i = 2*BodyTree.MaxBodiesNumber:-1:1 %Iterative over all the bodies discarding the fake ones
                if i <= l_B
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    if obj.BodiesInternal{i}.n ~= 0
                        k_b_1 = k_b - obj.BodiesInternal{i}.n + 1;
                        Dq(k_b_1:k_b) = obj.BodiesInternal{i}.D_;
                        %Prepare for the next iteration
                        k_b = k_b_1 - 1;
                    end
                end
            end
        end

        %Mass matrix.
        %NOTE: Current implementation has cost O(n^2)!
        function M = MassMatrix(obj, q)
            %Store value of gravity and set gravity to zero to compute the
            %mass matrix
            g_      = obj.g;
            obj.g   = zeros(3, 1, "like", g_);
            M       = zeros(obj.n, obj.n, "like", q);
            Zeron   = zeros(obj.n, 1, "like", q);
            Id      = eye(obj.n, "like", q);
            for i = 1:obj.n
                M(:, i) = obj.InverseDynamics(q, Zeron, Id(:, i));
            end
            %Symmetrize the result because of possible numerical
            %approximations
            M       = 1/2*(M + M');
            %Restrore gravity and compute the generalized forces
            obj.g   = g_;
        end

        %Gravity force
        function G = GravityForce(obj, q)
            Zeron = zeros(obj.n, 1, "like", q);
            G     = obj.InverseDynamics(q, Zeron, Zeron);
        end

        function C = ApparentForce(obj, q, dq)
            %Set gravity force to zero
            g_ = obj.g;
            obj.g = zeros(3, 1, "like", g_);
            %Compute the apparent forces
            C = obj.InverseDynamics(q, dq, zeros(obj.n, 1, "like", q));
            %Restore back gravity
            obj.g = g_;
        end
        
        
        %Computes the apparent (Coriolis/Centrifugal) matrix. We use the
        %approach of [Kawasaki et al., IEEE T-RA 1996]. The cost is
        %quadratic.
        function C = ApparentMatrix(obj, q, dq)
            %Set gravity force to zero
            g_      = obj.g;
            obj.g   = zeros(3, 1, "like", g_);
            %Preallocate the apparent matrix
            C       = zeros(obj.n, obj.n, "like", q);
            In      = eye(obj.n, "like", q);
            Zeron   = zeros(obj.n, 1, "like", q);
            %Compute the apparent forces
            for i = 1:obj.n
                C(:, i) = 1/2*( obj.InverseDynamics(q, dq + In(:, i), Zeron) ...
                               -obj.InverseDynamics(q, dq           , Zeron) ...
                               -obj.InverseDynamics(q, In(:, i)     , Zeron));
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
        

        function E = Energy(obj, q, dq, q_ref)
            % Compute the system energy.
            %q_ref is the reference configuration for the computation of
            %the elastic and gravitational energy. If it is not provided as
            %argument, 0 is assumed.
            switch nargin
                case 2
                    dq   = zeros(obj.n, 1, "like", q);
                    q_ref = zeros(obj.n, 1, "like", q);
                case 3
                    q_ref = zeros(obj.n, 1, "like", q);
            end
            
            %Kinetic energy
            E_kinetic = 1/2*dq'*obj.MassMatrix(q)*dq;
            %Potential energy
            E_elastic = 1/2*(q - q_ref)'*obj.K(q);
            if ~isnumeric(q)
                E_gravity = potential(obj.GravityForce(q), q, q_ref);
            else
                E_gravity = integral(@(s) (q - q_ref)'*obj.GravityForce(q_ref + s.*(q - q_ref)), 0, 1, 'ArrayValued', true);
            end
            %Overall energy
            E = E_kinetic + E_elastic + E_gravity;
        end

        function T = DirectKinematics(obj, q)
            %Compute the direct kinematics for each body. 
            % T is a cell array of transformation matrices from the base to the tip.
            T   = cell(obj.N_B, 1);
            T_i = obj.T0;
            obj = obj.TreeUpdate(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q));
            for i = 1:obj.N_B
                T_i  = T_i*obj.Joints{i}.T_*obj.Bodies{i}.T_;
                T{i} = T_i;
            end
        end
    end

    methods (Access = protected)
        
        %Implements the equlibrium equation
        function eq = EquilibriumEquation(obj, q, tau)
            eq = obj.GravityForce(q, 'double') + obj.K(q, 'double') - tau;
        end
        
    end
end

