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
        %Total number of actuators, required for code generation
        n_a  = 0;
        
        MassConditionNumber = 0;%TODO: Remove
    end
    
    properties (Access = protected)
        % Joints can be treated as bodies with no inertial parameters. Thus,
        % the augments the state by considering also the joints as bodies.
        BodiesInternal;
        N_B_Internal = 0;
        % Matrix that stores the start and end indexes of the configuration
        % vector for each body
        BodyConfigurationIndexes;
        % Matrix that stores the start and end indexes of the configuration
        % vector for each jointp
        JointConfigurationIndexes;
    end

    methods (Static)
        % Load a BodyTree object from its structure representation
        function obj = loadobj(S)
            % Create the cell arrays containing the joints and bodies of the tree
            Joints = cell(BodyTree.MaxBodiesNumber, 1);
            Bodies = cell(BodyTree.MaxBodiesNumber, 1);
            
            BodiesClasses = S.BodiesClasses;
            BodiesArray   = S.BodiesArray;
            JointsClasses = S.JointsClasses;
            JointsArray   = S.JointsArray;

            NB = S.N_B;
            for i = 1:BodyTree.MaxBodiesNumber
                if i <= NB
                    % Call the loadobj method of the body
                    loadobjBody     = str2func(string(BodiesClasses{i}) + ".loadobj");
                    Bodies{i}       = loadobjBody(BodiesArray{i});
                    % Call the loadobj method of the joint
                    loadobjJoint    = str2func(string(JointsClasses{i}) + ".loadobj");
                    Joints{i}       = loadobjJoint(JointsArray{i});
                else
                    Bodies{i} = 0;
                    Joints{i} = 0;
                end
            end

            % Build the class instance using the class name to automatically handle classes inerihthing from the BodyTree class that do not override this method
            ClassConstructor = str2func(S.ClassName);
            obj = ClassConstructor(Joints, Bodies);
            %obj = BodyTree(Joints, Bodies);

            % Restore the other parameters of the tree
            obj.T0 = S.T0;
            obj.g  = S.g;
            obj.MassConditionNumber = S.MassConditionNumber;
        end
    end

    methods
        % Store a structure representation of the bodytree
        function S = saveobj(obj)
            BodiesStructArray              = cellfun(@saveobj, obj.Bodies, "UniformOutput", false);
            BodiesClasses                  = cellfun(@class, obj.Bodies, "UniformOutput", false);
            JointsStructArray              = cellfun(@saveobj, obj.Joints, "UniformOutput", false);
            JointsClasses                  = cellfun(@class, obj.Joints, "UniformOutput", false);
            S = struct('n', obj.n, ...
                       'T0', obj.T0, ...
                       'g', obj.g, ...
                       'MassConditionNumber', obj.MassConditionNumber, ...
                       'N_B', obj.N_B, ...
                       'BodiesArray', {BodiesStructArray}, ...
                       'BodiesClasses', {BodiesClasses}, ...
                       'JointsArray', {JointsStructArray}, ...
                       'JointsClasses', {JointsClasses}, ...
                       'ClassName', class(obj), ...
                       'n_a', obj.n_a);% Number of actuators
        end
    end
    
    methods
        function obj = BodyTree(Joints, Bodies)
            %Construct a BodyTree consisting of joints and
            %bodies.
            %
            %Args:
                %    Joints ({:class:`Joint`}): Cell array of joints of the tree
                %    Bodies ({:class:`Body`}) : Cell array of bodies of the tree

            if ~iscolumn(Joints)
                Joints = Joints';
            end
            if ~iscolumn(Bodies)
                Bodies = Bodies';
            end
            
            % Compute the number of bodies and assign the joint and body
            % configuration indexes
            obj.N_B = 0;
            obj.n   = 0;
            obj.BodyConfigurationIndexes  = zeros(BodyTree.MaxBodiesNumber, 2);
            obj.JointConfigurationIndexes = zeros(BodyTree.MaxBodiesNumber, 2);
            prevIdx = 0;
            l_B     = length(Bodies);
            for i = 1:BodyTree.MaxBodiesNumber
                if i <= l_B
                    if ~isnumeric(Bodies{i})
                        obj.n = obj.n + Bodies{i}.n + Joints{i}.n;
                        obj.N_B = obj.N_B + 1;
                        % Assign the joint configuration indexes
                        if Joints{i}.n ~= 0
                            obj.JointConfigurationIndexes(i, 1:2) = [prevIdx + 1, prevIdx + Joints{i}.n];
                            % Store the index of the previous joint to use in the next iteration
                            prevIdx = obj.JointConfigurationIndexes(i, 2);
                        end
                        % Assign the body configuration indexes
                        if Bodies{i}.n ~= 0
                            obj.BodyConfigurationIndexes(i, 1:2)  = [prevIdx + 1, prevIdx + Bodies{i}.n];
                            % Store the index of the body to use in the next iteration
                            prevIdx = obj.BodyConfigurationIndexes(i, 2);
                        end
                    end
                end
            end
            % In a BodyTree class the number of actuators is the same as that of DOFs
            obj.n_a = obj.n;
            % Assign the joints and bodies
            % obj.Joints = cell(BodyTree.MaxBodiesNumber, 1);
            % obj.Bodies = cell(BodyTree.MaxBodiesNumber, 1);
            % if coder.target("MATLAB")% Augment the joints and bodies
            %      Joints = [Joints; cell(obj.MaxBodiesNumber-obj.N_B, 1)];
            %      Bodies = [Bodies; cell(obj.MaxBodiesNumber-obj.N_B, 1)];
            % end
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
       
        
        % Non parallelized version of the tree update
        function obj = TreeUpdate(obj, q, dq, ddq, options)
            %Update the state of the BodyTree. 
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of configuration variables
            
            % Arguments definition
            arguments
                obj (1, 1) BodyTree
                q {mustBeVector}
                dq {mustBeVector}
                ddq {mustBeVector}
                options.EvaluateKinematicTerms (1, 1) logical = true
                options.EvaluateInertialTerms (1, 1) logical  = true
                options.EvaluateExternalForces (1, 1) logical = true
            end
            
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
            l_B = length(obj.Bodies);
            for i = 1:BodyTree.MaxBodiesNumber
                if i <= l_B
                    if isnumeric(obj.Bodies{i}) || isnumeric(obj.Joints{i})
                        continue;
                    end
                    % JOINT UPDATE
                    % Get the DOFs associated with the joint
                    n_j   = obj.Joints{i}.n;
                    if n_j ~= 0
                        % Get the indexes associated with the joint
                        q_start = obj.JointConfigurationIndexes(i, 1);
                        q_end   = obj.JointConfigurationIndexes(i, 2);
                        % Update the joint data
                        obj.Joints{i}.Update(q(q_start:q_end), dq(q_start:q_end), ddq(q_start:q_end), "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                                                                                      "EvaluateInertialTerms" ,  options.EvaluateInertialTerms, ...
                                                                                                      "EvaluateExternalForces", options.EvaluateExternalForces);
                    %else
                    %    obj.Joints{i}.Update([], [], [], "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                    %                                     "EvaluateInertialTerms" , options.EvaluateInertialTerms, ...
                    %                                     "EvaluateExternalForces", options.EvaluateExternalForces);
                    end
                    
                    % BODY UPDATE
                    % Get the DOFs associated with the body
                    n_b   = obj.Bodies{i}.n;
                    if n_b ~= 0
                        % Get the indexes associated with the body
                        q_start = obj.BodyConfigurationIndexes(i, 1);
                        q_end   = obj.BodyConfigurationIndexes(i, 2);
                        % Update the body data
                        obj.Bodies{i}.Update(q(q_start:q_end), dq(q_start:q_end), ddq(q_start:q_end), "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                                                                                      "EvaluateInertialTerms" ,  options.EvaluateInertialTerms, ...
                                                                                                      "EvaluateExternalForces", options.EvaluateExternalForces);
                    %else
                    %    obj.Bodies{i}.Update([], [], [], "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                    %                                     "EvaluateInertialTerms" , options.EvaluateInertialTerms, ...
                    %                                     "EvaluateExternalForces", options.EvaluateExternalForces);
                    end
                end
            end
        end


        % Parallilzed verion of the tree update
        function obj = ParTreeUpdate(obj, q, dq, ddq, options)
            %Update the state of the BodyTree. 
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of configuration variables
            
            % Arguments definition
            arguments
                obj (1, 1) BodyTree
                q {mustBeVector}
                dq {mustBeVector}
                ddq {mustBeVector}
                options.EvaluateKinematicTerms (1, 1) logical = true
                options.EvaluateInertialTerms (1, 1) logical  = true
                options.EvaluateExternalForces (1, 1) logical = true
            end
            
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
            l_B  = length(obj.Bodies);
            B    = obj.Bodies;
            J    = obj.Joints;
            JIdx1 = obj.JointConfigurationIndexes(:, 1);
            JIdx2 = obj.JointConfigurationIndexes(:, 2);
            BIdx1 = obj.BodyConfigurationIndexes(:, 1);
            BIdx2 = obj.BodyConfigurationIndexes(:, 2);
            parfor i = 1:BodyTree.MaxBodiesNumber
                Ji    = J{i};
                Bi    = B{i};
                    
                if i <= l_B
                    if isnumeric(B{i}) || isnumeric(J{i})
                        continue;
                    end
                    % JOINT UPDATE
                    % Get the DOFs associated with the joint
                    n_j   = Ji.n;
                    if n_j ~= 0
                        % Get the indexes associated with the joint
                        q_start = JIdx1(i);
                        q_end   = JIdx2(i);
                        % Update the joint data
                        Ji.Update(q(q_start:q_end), dq(q_start:q_end), ddq(q_start:q_end), "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                                                                                      "EvaluateInertialTerms" ,  options.EvaluateInertialTerms, ...
                                                                                                      "EvaluateExternalForces", options.EvaluateExternalForces);
                    else
                        Ji.Update([], [], [], "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                                         "EvaluateInertialTerms" ,  options.EvaluateInertialTerms, ...
                                                         "EvaluateExternalForces", options.EvaluateExternalForces);
                    end
                    
                    % BODY UPDATE
                    % Get the DOFs associated with the body
                    n_b   = Bi.n;
                    if n_b ~= 0
                        % Get the indexes associated with the body
                        q_start = BIdx1(i);
                        q_end   = BIdx2(i);
                        % Update the body data
                        Bi.Update(q(q_start:q_end), dq(q_start:q_end), ddq(q_start:q_end), "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                                                                                      "EvaluateInertialTerms" ,  options.EvaluateInertialTerms, ...
                                                                                                      "EvaluateExternalForces", options.EvaluateExternalForces);
                    else
                        Bi.Update([], [], [], "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                                         "EvaluateInertialTerms" ,  options.EvaluateInertialTerms, ...
                                                         "EvaluateExternalForces", options.EvaluateExternalForces);
                    end

                end

                B{i} = Bi;
                J{i} = Ji;
            end
            for i = 1:BodyTree.MaxBodiesNumber
                obj.Bodies{i} = B{i};
                obj.Joints{i} = J{i};
                if mod(i, 2) == 0
                    obj.BodiesInternal{i} = B{i};
                else
                    obj.BodiesInternal{i} = J{i};
                end
            end
        end

      
        function tau = InverseDynamics(obj, q, dq, ddq, options)
            % Evaluate the inverse dynamics using the Generalized Inverse Dynamics (GID) algorithm.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of configuration variables

            arguments (Input)
                obj (1, 1) BodyTree
                q   (:, 1)
                dq  (:, 1)
                ddq (:, 1)
                options.EvaluateKinematicTerms (1, 1) logical  = true
                options.EvaluateInertialTerms  (1, 1) logical  = true
                options.EvaluateExternalForces (1, 1) logical  = true
            end
            
            % Update the tree
            obj = obj.TreeUpdate(q, dq, ddq, ...
                                    "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                    "EvaluateInertialTerms",  options.EvaluateInertialTerms, ...
                                    "EvaluateExternalForces", options.EvaluateExternalForces);
            % Run the GID algorithm,  q is only passed to
            % retreive its type.
            tau = obj.Kane_aux(obj.g, q);
            
            % Add the other external forces, i.e., damping and elastic force, if requested.
            if options.EvaluateExternalForces == true
                tau = tau ...
                    + obj.K(q, "EvaluateKinematicTerms", false, "EvaluateInertialTerms", false) ...
                    + obj.D(q, dq, "EvaluateKinematicTerms", false, "EvaluateInertialTerms", false);
            end
            
        end

        function [Mq, tau] = UnifiedInverseDynamics(obj, q, dq, ddq, options)
            % Evaluate the unified inverse dynamics using the Generalized Inverse Dynamics (GID) algorithm.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of configuration variables

            arguments (Input)
                obj (1, 1) BodyTree
                q   (:, 1)
                dq  (:, 1)
                ddq (:, 1)
                options.EvaluateKinematicTerms (1, 1) logical  = true
                options.EvaluateInertialTerms  (1, 1) logical  = true
                options.EvaluateExternalForces (1, 1) logical  = true
            end
            
            % Update the tree
            obj = obj.TreeUpdate(q, dq, ddq, ...
                                    "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                    "EvaluateInertialTerms",  options.EvaluateInertialTerms, ...
                                    "EvaluateExternalForces", options.EvaluateExternalForces);
            % Run the GID algorithm,  q is only passed to
            % retreive its type.
            [Mq, tau] = obj.UnifiedKane(obj.g, q);

            Mq = 1/2*(Mq + Mq');
            
            % Add the other external forces, i.e., damping and elastic force, if requested.
            if options.EvaluateExternalForces == true
                tau = tau ...
                    + obj.K(q, "EvaluateKinematicTerms", false, "EvaluateInertialTerms", false) ...
                    + obj.D(q, dq, "EvaluateKinematicTerms", false, "EvaluateInertialTerms", false);
            end
            
        end

        
        function [Mq, tau] = ParUnifiedInverseDynamics(obj, q, dq, ddq, options)
            % Evaluate the unified inverse dynamics using the Generalized Inverse Dynamics (GID) algorithm.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of configuration variables

            arguments (Input)
                obj (1, 1) BodyTree
                q   (:, 1)
                dq  (:, 1)
                ddq (:, 1)
                options.EvaluateKinematicTerms (1, 1) logical  = true
                options.EvaluateInertialTerms  (1, 1) logical  = true
                options.EvaluateExternalForces (1, 1) logical  = true
            end
            
            % Update the tree
            obj = obj.ParTreeUpdate(q, dq, ddq, ...
                                    "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                    "EvaluateInertialTerms",  options.EvaluateInertialTerms, ...
                                    "EvaluateExternalForces", options.EvaluateExternalForces);
            % Run the GID algorithm,  q is only passed to
            % retreive its type.
            [Mq, tau] = obj.UnifiedKane(obj.g, q);

            Mq = 1/2*(Mq + Mq');
            
            % Add the other external forces, i.e., damping and elastic force, if requested.
            if options.EvaluateExternalForces == true
                tau = tau ...
                    + obj.K(q, "EvaluateKinematicTerms", false, "EvaluateInertialTerms", false) ...
                    + obj.D(q, dq, "EvaluateKinematicTerms", false, "EvaluateInertialTerms", false);
            end
            
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
            G  = obj.GravityForce(q);
            C  = obj.ApparentForce(q, dq);% We evaluate last C so that both elastic and damping forces can be computed withou running again all the computations
            % Overall forces
            f  =  -C - G - obj.K(q, "EvaluateKinematicTerms", false, "EvaluateInertialTerms", false) - obj.D(q, dq, "EvaluateKinematicTerms", false, "EvaluateInertialTerms", false) + tau;
            % Perform scaling to improve simulation accuracy
            % Compute the acceleration
            ddq = M\f;
        end

        function Kq = K(obj, q, options)
            %Evaluate the generalized elastic force.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
            
            % Define arguments
            arguments
                obj                (1, 1) BodyTree
                q                  (:, 1)
                options.EvaluateKinematicTerms (1, 1) logical = true
                options.EvaluateInertialTerms  (1, 1) logical = true
            end
            % Update the tree, but only the stiffness force 
            if nargin == 2
                obj.TreeUpdate(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q), "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                                                                          "EvaluateExternalForces", true, ...
                                                                                          "EvaluateInertialTerms", options.EvaluateInertialTerms);
            end
            % switch nargin
            %     case 1
            %         q = zeros(obj.n, 1);
            %     case 2
            %         obj.TreeUpdate(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q));
            % end
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

        function Dq = D(obj, q, dq, options)
            %Evaluate the generalized damping force.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of configuration variables
            
            % Define arguments
            arguments
                obj (1, 1) BodyTree
                q   (:, 1)
                dq  (:, 1)
                options.EvaluateKinematicTerms (1, 1) logical = true
                options.EvaluateInertialTerms  (1, 1) logical = true
            end

            % Update the tree only if q and dq are passed as arguments
            if nargin >= 2
                obj.TreeUpdate(q, dq, zeros(obj.n, 1, "like", q), "EvaluateKinematicTerms", options.EvaluateKinematicTerms, ...
                                                                  "EvaluateExternalForces", true, ...
                                                                  "EvaluateInertialTerms", options.EvaluateInertialTerms);
            end
            % switch nargin
            %     case 1
            %         q = zeros(obj.n, 1);
            %     case 2
            %         obj.TreeUpdate(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q));
            %     case 3
            %         obj.TreeUpdate(q, dq, zeros(obj.n, 1, "like", q));
            % end
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
                M(:, i) = obj.InverseDynamics(q, Zeron, Id(:, i), "EvaluateKinematicTerms", true, "EvaluateInertialTerms", true, "EvaluateExternalForces", false);
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
            G     = obj.InverseDynamics(q, Zeron, Zeron, "EvaluateKinematicTerms", true, "EvaluateInertialTerms", true, "EvaluateExternalForces", false);
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
            C = obj.InverseDynamics(q, dq, zeros(obj.n, 1, "like", q), "EvaluateKinematicTerms", true, "EvaluateInertialTerms", true, "EvaluateExternalForces", false);
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
        
        function [q_eq, f_val] = EquilibriumConfiguration(obj, q0, tau, options)
            %Find an equilibrium configuration by solving numerically the equilibrium equations.
            %
            %Args:
                %    q0  ([double]): Initial guess for the equlibrium
                %    tau ([double]): Generalized actuation force
            arguments
                obj (1, 1) BodyTree
                q0  (:, 1) double
                tau (:, 1) double
                options.FunctionTolerance   (1, 1) double = 1e-6
                options.StepTolerance       (1, 1) double = 1e-6
            end
            opt_options = optimoptions('fsolve', ...
                                       'FunctionTolerance', options.FunctionTolerance, ...
                                       'StepTolerance', options.StepTolerance);
            [q_eq, f_val] = fsolve(@(q) obj.EquilibriumEquation(q, tau), q0, opt_options);
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
            %   idx           ([int]): Ordered array of body indexes indicating the bodies for which the direct kinematics has to be evaluated
            %Return:
            %   ([double], [sym]): Homogeneous transformation matrices for the body specified by idx.   

            % Arguments definition
            arguments (Input)
                obj (1, 1) BodyTree
                q {mustBeVector}    = zeros(obj.n, 1)
                idx {mustBeVector}  = linspace(1, obj.N_B, obj.N_B);
            end

            arguments (Output)
                %T (:, 4)
                T (4, 4, :)
            end

            % Update the state of the kinematic tree
            Zeron   = zeros(obj.n, 1, "like", q);
            obj     = obj.TreeUpdate(q, Zeron, Zeron, "EvaluateKinematicTerms", true, ...
                                                      "EvaluateInertialTerms", false, ...
                                                      "EvaluateExternalForces", false);

            % Length of the index vector
            idxLength = length(idx);

            % Modify the index to account for the fact that internally the joints are modeled as bodies
            idx = 2*idx;

            % Check that the input vectors are in column format.
            if ~iscolumn(q)
                q = q';
            end

            % T is matrix of vertically stacked 4x4 transformation matrices
            %T   = repmat(eye(4, 'like', q), length(idx), 1);
            T   = repmat(eye(4, 'like', q), 1, 1, length(idx));
            
            % Variables initialization
            T_i         = obj.T0;
            lastIdx     = idx(end);
            N_B_        = obj.N_B_Internal;
            % Compute the direct kinematics
            j           = 1;
            for i = 1:2*BodyTree.MaxBodiesNumber % Iterative over all the augmented bodies
                if i <= N_B_
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end

                    % Check if the index is zero, i.e., the base.
                    if idx(j) == 0
                        %T(1+4*(j-1):4*j, 1:4)   = T_i;
                        T(1:4, 1:4, j)          = T_i;
                        j                       = j + 1;
                        if j > idxLength % Verify if idx = [0]; in case exit to avoid access out of range
                            break;
                        end
                    end
                    
                    % Compute the transformation matrix from the base to
                    % the current body
                    T_i  = T_i*obj.BodiesInternal{i}.T_;
                    
                    % Check if the current body is in the index list and in
                    % case assign T_i as output
                    if i == idx(j)
                        %T(1+4*(j-1):4*j, 1:4)   = T_i;
                        T(1:4, 1:4, j)          = T_i;
                        j                       = j + 1;
                    end
    
                    % If we have reached the last body, break the loop
                    if i == lastIdx
                        break;
                    end
                end
            end
        end

        function [q, converged, e] = InverseKinematics(obj, T, options)
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

            arguments
                obj 
                T
                options.BodyIndexes             = 1:obj.N_B
                options.InitialGuess            = zeros(obj.n, 1)
                options.TaskFlags               = ones(obj.N_B*6, 1)
                options.MaxIterationNumber      = 4
                options.AngularErrorThreshold   = 1e-3
                options.LinearErrorThreshold    = 1e-2
                options.UseGradientDescent      = false
                options.GradientDescentStepSize = 1
                options.ErrorWeight             = 1
            end
            
            % Variable assignment
            idx         = options.BodyIndexes;
            q0          = options.InitialGuess;
            task_flags  = options.TaskFlags;
            N           = options.MaxIterationNumber;
            AngularErrorThsd = options.AngularErrorThreshold;
            LinearErrorThsd  = options.LinearErrorThreshold;
            UseGradientDescent = options.UseGradientDescent;
            GradientDescentStepSize = options.GradientDescentStepSize;
            ErrorWeight = options.ErrorWeight;

            % Store useful variables
            idxLength = length(idx);
            taskFlagsLength = length(task_flags);

            % Convert the task flags to an index vector
            task_flags_l = logical(task_flags == 0);
            
            % Check that the dimensions of T and idx are consistent. 
            if floor(size(T, 1)/4) ~= idxLength
                error("The size of T and idx is not consistent.");
            end
            % Check that the dimensions of idx and task_flags are consistent
            if 6*idxLength ~= taskFlagsLength
                error("The length of the body indexes is not consitent with the length of the task flags.");
            end
            % Check that the dimension of the weight matrix is consistent
            % with the dimension of the task vector
            [ew1, ew2] = size(ErrorWeight);
            if ew1 ~= ew2
                error("The weight matrix must be square.");
            end
            if ew1 ~= 1
                if taskFlagsLength ~= ew1
                    error("The size of the weight matrix must be equal the length of the task flags.");
                end
            end

            % Preallocate the output for code generation
            q         = zeros(obj.n, 1, "like", q0);% DO NOT REMOVE THIS!
            q         = q0;
            converged = 0;

            % Allocate iteration variables
            e         = Inf*ones(6*idxLength, 1);
            e_thsd    = ones(idxLength, 1);%If contains all zeros, the configuration satisfies all the constraints

            ErrorWeight12 = mpower(ErrorWeight, 0.5);
            
            % Run the Newton algorithm as given in Linch and Park, Modern Robotics
            for i = 1:N
                % Evaluate the direct kinematics in the current configuration
                T_q         = obj.DirectKinematics(q, idx);
                % Evaluate the body Jacobian in the current configuration and reshape it as vertically stacked jacobians
                J_q         = reshape(permute(obj.BodyJacobian(q, idx), [1, 3, 2]), 6*idxLength, []);
                
                % Iterate over all the requried bodies
                for j = 1:idxLength
                    % Evaluate the desired configuration in the current frame
                    %T_qd_j = invTransformation(T_q(1+4*(j-1):4*j, 1:4))*T(1+4*(j-1):4*j, 1:4);
                    T_qd_j = invTransformation(T_q(1:4, 1:4, j))*T(1+4*(j-1):4*j, 1:4);
                    
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
                    % Update the configuration using Gradient Descent
                    if UseGradientDescent == true
                        q = q + GradientDescentStepSize*(J_q')*ErrorWeight*e;
                    else% Update the configuration using weighted Newton-Rapson
                        %q = q + J_q\e;
                        % Use weighted pseudoinverse
                        q = q + (ErrorWeight12*J_q)\(ErrorWeight12*e);
                    end
                end
            end

            % Display a warning if the convergence was not achieved
            %if ~converged
            %    warning("Convergence not achived. The error norm is: " + sprintf("%.4f", norm(e)));
            %end

        end

        function plot(obj, q, options)
            % Plot the robot in the current configuration
            arguments
                obj                           (1, 1) BodyTree
                q                             (:, 1) double = zeros(obj.n, 1)
                options.Color                 (:, 3) double = [0 160 219]./256;
                options.FaceAlpha             (:, 1) double = 1
                options.LineStyle             (:, 1) string = "none";% Possible values: "-" | "--" | ":" | "-." | "none"
                options.TransformationMatrix  (4, 4) double = zeros(4, 4)
            end
            % Get the size of the plotting properties
            ColorLength     = size(options.Color, 1);
            FaceAlphaLength = size(options.FaceAlpha, 1);
            LineStyleLength = size(options.LineStyle, 1);
            
            % Utility function to circullary wrap vectors indexes during iteration
            wrapN = @(i, N) (1 + mod(i-1, N));
            

            % Update the status of the tree
            obj.TreeUpdate(q, zeros(obj.n, 1), zeros(obj.n, 1));

            % Iterate over the bodies, each body must define a plot method that handles its plots for their configuration
            if all(all(options.TransformationMatrix == 0))
                T0 = obj.T0;
            else
                T0 = options.TransformationMatrix;
            end
            
            % Counter for the properties
            j     = 1;
            for i = 1:obj.N_B
                % Update the transformation matrix using the joint data
                T0 = T0*obj.Joints{i}.T_;
                % Get also the transform induced by the body, since some classes might alter their state during plot we save it
                TB = obj.Bodies{i}.T_;
                % Plot the body
                n_b   = obj.Bodies{i}.n;
                % Check if the current bodies implements the plot method
                if ismethod(obj.Bodies{i}, "plot")
                    if n_b ~= 0
                        % Get the indexes associated with the body
                        q_start = obj.BodyConfigurationIndexes(i, 1);
                        q_end   = obj.BodyConfigurationIndexes(i, 2);
                        % Update the body data
                        obj.Bodies{i}.plot(q(q_start:q_end), "BaseTransformation", T0, ...
                                                             "Color", options.Color(wrapN(j, ColorLength), 1:3), ...
                                                             "FaceAlpha", options.FaceAlpha(wrapN(j, FaceAlphaLength)), ...
                                                             "LineStyle", options.LineStyle(wrapN(j, LineStyleLength)));
                    else
                        obj.Bodies{i}.plot("BaseTransformation", T0, ...
                                           "Color", options.Color(wrapN(j, ColorLength), 1:3), ...
                                           "FaceAlpha", options.FaceAlpha(wrapN(j, FaceAlphaLength)), ...
                                           "LineStyle", options.LineStyle(wrapN(j, LineStyleLength)));
                    end

                    % Update the counter for the properties
                    j = j + 1;
                end
                % Update the transformation matrix by accounting for the transform of the body
                T0 = T0*TB;
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
            
            % Arguments definition
            arguments (Input)
                obj (1, 1) BodyTree
                q {mustBeVector}    = zeros(obj.n, 1)
                idx {mustBeVector}  = linspace(1, obj.N_B, obj.N_B)
            end

            arguments (Output)
                J (6, :, :)
            end

            % Update the BodyTree
            Zeron   = zeros(obj.n, 1, 'like', q);
            obj     = obj.TreeUpdate(q, Zeron, Zeron, "EvaluateKinematicTerms", true, ...
                                                      "EvaluateExternalForces", false, ...
                                                      "EvaluateInertialTerms", false);

            % Modify the index to account for the fact that internally the joints are modeled as bodies
            idx = 2*idx;

            % Check that the input vectors are in column format.
            if ~iscolumn(q)
                q = q';
            end

            % Define auxiliary variables
            % The joints are treated as massless bodies thus we augment the
            % body dimensions
            N_B_ = obj.N_B_Internal;
            
            % Define the output
            %J      = zeros(length(idx)*6, obj.n, 'like', q);
            J      = zeros(6, obj.n, length(idx), 'like', q);
            
            
            % Auxiliary variables for the iteration
            J_i        = zeros(6, obj.n, 'like', q);
            q_idx      = 1;
            lastIdx    = idx(end);
            
            % Check if the first index is the one for the base, i.e., is 0
            if idx(1) == 0
                j          = 2;
            else
                j          = 1;
            end
            % Check if j exceeds the length of idx, i.e., idx = [0], in
            % case return.
            if j > length(idx)
                return;
            end
            
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

                    % Update the Jacobian for the previous variables by rotating it in the
                    % frame of the current body
                    nBody                       = obj.BodiesInternal{i}.n;
                    if q_idx ~= 1
                        J_i(4:6, 1:q_idx-1)       = real(R_i_T*(J_i(4:6, 1:q_idx-1) - skew(t_i)*J_i(1:3, 1:q_idx-1)));
                        J_i(1:3, 1:q_idx-1)       = real(R_i_T*J_i(1:3, 1:q_idx-1));
                    end
                    J_i(1:6, q_idx:q_idx+nBody-1) = [omega_par_i; v_par_i];

                    % Save the value if i hits the current body value
                    if i == idx(j)
                        %J(1 + 6*(j-1):6*j, 1:obj.n) = J_i;
                        J(1:6, 1:obj.n, j) = J_i;
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
            
        end

        function A = ActuationMatrix(obj, q)
        % Method to be overridden by custom classes to implement the actuation matrix computation.
            arguments (Input)
                obj (1, 1) BodyTree
                q   (:, 1)
            end

            arguments (Output)
                A (:, :)
            end
            
            % By default the actuation matrix is the identity
            A = eye(obj.n);
        end

        % Getter methods.
        
        
        function Bidx = getBodyConfigurationIndex(obj, idx)
            % Get the indexes of the configuration vector for the bodies
            % specified by idx.
            %
            %Args:
            %   idx   ([double, int]): Bodies for which the configuration
            %   indexes have to be computed.
            %Return:
            %   {[double], [sym]}: length(idx) x 2 matrix where each row
            %   specifies the start and end index of the corresponding body
            %   in the idx vector. If idx is not specified, then the
            %   function returns the indexes for all the bodies in the
            %   tree.
            
            % If the method is called without indexes, return the entire
            % vector of indexes
            if nargin == 1
                Bidx = obj.BodyConfigurationIndexes;
                return;
            end
            Bidx = obj.BodyConfigurationIndexes(idx, 1:2);
        end

        function Jidx = getJointConfigurationIndex(obj, idx)
            % Get the indexes of the configuration vector for the joints
            % specified by idx.
            %
            %Args:
            %   idx   ([double, int]): Joints for which the configuration
            %   indexes have to be computed.
            %Return:
            %   {[double], [sym]}: length(idx) x 2 matrix where each row
            %   specifies the start and end index of the corresponding
            %   joint in the idx vector. If idx is not specified, then the
            %   function returns the indexes for all the joints in the
            %   tree.
            if nargin == 1
                Jidx = obj.JointConfigurationIndexes;
                return;
            end
            Jidx = obj.JointConfigurationIndexes(idx, 1:2);
        end
    end

    methods (Access = protected)
        
        function eq = EquilibriumEquation(obj, q, tau)
            %Evaluate the equilibrium equation.
            %
            %Args:
                %    q      ([double], [sym]): Configuration variables
                %    tau    ([double], [sym]): Generalized actuation force
            eq = obj.GravityForce(q) + obj.K(q) - obj.ActuationMatrix(q)*tau;
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
    
    
    
        function [Mq, tau] = UnifiedKane(obj, g, q_type)
            %Compute the unified generalized inverse dynamics using the recursion
            %formulas.
            %Args:
                %    g      ([double], [sym]): Gravitational force vector expressed in the base frame
                %    q_type ([double], [sym]): Dummy variable used to retrieve the type of the configuration variables
            
            % Define auxiliary variables
            % The joints are treated as massless bodies thus we augment the
            % body dimensions
            N_B_            = obj.N_B_Internal;
            
            % Define the output
            Mq              = zeros(obj.n, obj.n, "like", q_type);
            tau             = zeros(obj.n, 1, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the linear velocity of Body i ( origin of frame {S_i} )
            v               = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the angular velocity of Body i
            omega           = zeros(3, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the linear acceleration of Body i ( origin of frame {S_i} )
            a               = zeros(3, N_B_, "like", q_type);
            % 3 x n x N_B tensor that stores the Jacobian of the linear acceleration of Body i w.r.t. \ddot{q}
            Ja              = zeros(3, obj.n, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the angular acceleration of Body i
            domega          = zeros(3, N_B_, "like", q_type);
            % 3 x n x N_B tensor that stores the Jacobian of the angular acceleration of Body i w.r.t. \ddot{q}
            Jdomega         = zeros(3, obj.n, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores the linear acceleration of CoM_i
            a_com           = zeros(3, N_B_, "like", q_type);
            % 3 x n x N_B tensor that stores the Jacobian of the CoM acceleration of Body i w.r.t. \ddot{q}
            Ja_com          = zeros(3, obj.n, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores vector Gamma_i
            Gamma           = zeros(3, N_B_, "like", q_type);
            % 3 x n x N_B tensor that stores the Jacobian of Gamma_i w.r.t. \ddot{q}
            JGamma          = zeros(3, obj.n, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores vector Omega_i
            Omega           = zeros(3, N_B_, "like", q_type);
            % 3 x n x N_B tensor that stores the Jacobian of Omega_i w.r.t. \ddot{q}
            JOmega          = zeros(3, obj.n, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores vector M_i
            M               = zeros(3, N_B_, "like", q_type);
            % 3 x n x N_B tensor that stores the Jacobian of M_i w.r.t. \ddot{q}
            JM              = zeros(3, obj.n, N_B_, "like", q_type);
            % 3 x N_B matrix whose i-th column stores vector N_i
            N               = zeros(3, N_B_, "like", q_type);
            % 3 x n x N_B tensor that stores the Jacobian of N_i w.r.t. \ddot{q}
            JN              = zeros(3, obj.n, N_B_, "like", q_type);
            % Auxiliary vector to represent the velocity of the preceding body
            v_i_1           = zeros(3, 1, "like", q_type);
            omega_i_1       = zeros(3, 1, "like", q_type);
            % Auxiliary vectors to represent the accelerations of the preceding body
            a_i_1           = -g;
            domega_i_1      = zeros(3, 1, "like", q_type);
            M_i_1           = zeros(3, 1, "like", q_type);
            N_i_1           = zeros(3, 1, "like", q_type);
            % Auxiliary vectors to represent the Jacobians of the accelerations of the preceding body
            Ja_i_1          = zeros(3, obj.n, "like", q_type);
            Jdomega_i_1     = zeros(3, obj.n, "like", q_type);
            JM_i_1          = zeros(3, obj.n, "like", q_type);
            JN_i_1          = zeros(3, obj.n, "like", q_type);
            % Auxiliary variable that counts the number of DOFs already processed
            n_i             = 0;
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
                                          -skew(t_i)*domega_i_1 + ...
                                          skew(omega_i_1)*(skew(omega_i_1)*t_i + v_rel_i) + ...
                                          skew(omega_i_1)*v_rel_i + ...
                                          dv_rel_i));
                    domega(:, i) = real(R_i_T*(domega_i_1 + skew(omega_i_1)*omega_rel_i + domega_rel_i));
                    a_com(:, i)  = real(a(:, i) - skew(p_com_i)*domega(:, i) + cross(omega(:, i), cross(omega(:, i), p_com_i) + v_com_rel_i) +...
                                    cross(omega(:, i), v_com_rel_i) + a_com_rel_i);
                    
                    n_Body       = obj.BodiesInternal{i}.n;
                    Ja(:, :, i)  = real(R_i_T*(Ja_i_1 -skew(t_i)*Jdomega_i_1) ...
                                       + [zeros(3, n_i), obj.BodiesInternal{i}.v_par_, zeros(3, obj.n - n_i - n_Body)]);
                    Jdomega(:, :, i) = real(R_i_T*(Jdomega_i_1) +  [zeros(3, n_i), obj.BodiesInternal{i}.omega_par_, zeros(3, obj.n - n_i - n_Body)]);
                    Ja_com(:, :, i)   = real(Ja(:, :, i) - skew(p_com_i)*Jdomega(:, :, i) + [zeros(3, n_i), obj.BodiesInternal{i}.grad_v_com_', zeros(3, obj.n - n_i - n_Body)]);
                    
                    
                    
                    % Step 2
                    Gamma(:, i) =real(a_com(:, i)*obj.BodiesInternal{i}.m_ +...
                                    2*cross(omega(:, i), obj.BodiesInternal{i}.int_dr_) +...
                                    obj.BodiesInternal{i}.int_ddr_);
                    Omega(:, i) = real(obj.BodiesInternal{i}.I_*domega(:, i) + cross(omega(:, i), obj.BodiesInternal{i}.I_*omega(:, i)) +...
                                    obj.BodiesInternal{i}.J_*omega(:, i)+...
                                    cross(omega(:, i), obj.BodiesInternal{i}.int_r_X_dr_) + ...
                                    obj.BodiesInternal{i}.int_r_X_ddr_);

                    JGamma(:, :, i) = real(Ja_com(:, :, i)*obj.BodiesInternal{i}.m_ + [zeros(3, n_i), obj.BodiesInternal{i}.Jint_ddr_, zeros(3, obj.n - n_i - n_Body)]);
                    JOmega(:, :, i) = real(obj.BodiesInternal{i}.I_*Jdomega(:, :, i) + [zeros(3, n_i), obj.BodiesInternal{i}.Jint_r_X_ddr_, zeros(3, obj.n - n_i - n_Body)]);
                    
                    
                    % Update the iteration variables
                    v_i_1 = v(:, i);
                    omega_i_1 = omega(:, i);
                    a_i_1 = a(:, i);
                    domega_i_1 = domega(:, i);

                    Ja_i_1      = Ja(:, :, i);
                    Jdomega_i_1 = Jdomega(:, :, i);
                    n_i         = n_i + n_Body;
                end
            end
            
            % ********************************************************
            % ********************* Backward step ********************
            % ********************************************************
            % Index for tau vector
            idx_tau = obj.n;
            n_i     = obj.n;
            for i = 2*BodyTree.MaxBodiesNumber:-1:1
                if i <= N_B_
                    if isnumeric(obj.BodiesInternal{i})
                        continue;
                    end
                    % Step 3
                    if i == N_B_
                        M(:, i) = real(Gamma(:, i));
                        N(:, i) = real(Omega(:, i) + skew(obj.BodiesInternal{i}.p_com_)*Gamma(:, i));
                        
                        JM(:, :, i) = real(JGamma(:, :, i));
                        JN(:, :, i) = real(JOmega(:, :, i) + skew(obj.BodiesInternal{i}.p_com_)*JGamma(:, :, i));
                    else
                        if ~isnumeric(obj.BodiesInternal{i+1})
                            R_i_1   = real(obj.BodiesInternal{i+1}.T_(1:3, 1:3));
                            t_i_1   = real(obj.BodiesInternal{i+1}.T_(1:3, 4));
                            M(:, i) = real(Gamma(:, i) + R_i_1*M_i_1);
                            N(:, i) = real(Omega(:, i) + skew(obj.BodiesInternal{i}.p_com_)*Gamma(:, i) + R_i_1*N_i_1 + skew(t_i_1)*(R_i_1*M_i_1));

                            JM(:, :, i) = real(JGamma(:, :, i) + R_i_1*JM_i_1);
                            JN(:, :, i) = real(JOmega(:, :, i) + skew(obj.BodiesInternal{i}.p_com_)*JGamma(:, :, i) + R_i_1*JN_i_1 + skew(t_i_1)*(R_i_1*JM_i_1));
                        end
                    end
                    
                    % Step 4
                    if obj.BodiesInternal{i}.n ~= 0
                        n_Body = obj.BodiesInternal{i}.n;

                        % Update the iteration variable for the DOFs
                        n_i         = n_i - n_Body;

                        pi_star                         = real(obj.BodiesInternal{i}.grad_int_dr_*a_com(:, i) +...
                                                                     obj.BodiesInternal{i}.grad_int_r_X_dr_*domega(:, i) +...
                                                                     arrayfun(@(j) -(1/2)*omega(:, i)'*obj.BodiesInternal{i}.grad_J_(:, :, j)*omega(:, i), (1:size(obj.BodiesInternal{i}.grad_J_, 3))')+...
                                                                     2*obj.BodiesInternal{i}.int_dr_X_pv_r_*omega(:, i) +...
                                                                     obj.BodiesInternal{i}.int_pv_r_O_dd_r_ +...
                                                                     obj.BodiesInternal{i}.grad_v_com_*Gamma(:, i));
                        J_pi_star                       = real(obj.BodiesInternal{i}.grad_int_dr_*Ja_com(:, :, i) +...
                                                                     obj.BodiesInternal{i}.grad_int_r_X_dr_*Jdomega(:, :, i) +...
                                                                     [zeros(n_Body, n_i), obj.BodiesInternal{i}.Jint_pv_r_O_dd_r_, zeros(n_Body, obj.n - n_i - n_Body)] +...
                                                                     obj.BodiesInternal{i}.grad_v_com_*JGamma(:, :, i));
                    
                        tau(idx_tau - n_Body + 1:idx_tau) = real(pi_star +...
                                                                 (obj.BodiesInternal{i}.v_par_')*M(:, i) +...
                                                                 (obj.BodiesInternal{i}.omega_par_')*N(:, i));
                    
                        Mq(idx_tau - n_Body + 1:idx_tau, 1:obj.n) = real( J_pi_star +...
                                                                     (obj.BodiesInternal{i}.v_par_')*JM(:, :, i) +...
                                                                     (obj.BodiesInternal{i}.omega_par_')*JN(:, :, i));

                    end
                    % Update the iteration variables
                    M_i_1 = M(:, i);
                    N_i_1 = N(:, i);

                    JM_i_1 = JM(:, :, i);
                    JN_i_1 = JN(:, :, i);
                    idx_tau = idx_tau - obj.BodiesInternal{i}.n;
                end
            end
        end
    
    end
end

