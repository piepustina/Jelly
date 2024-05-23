classdef LVPBodyTreeJAct < BodyTree
    %LVPBODYTREEJACT Test class used for defining custom actuation matrix
    
    properties (Access = private)
        ActuatorPosition  (:, 3)
        ActuatorBody      (:, 1)
        ActuatorDirection (:, 3)
        ClosestBodyCentroid (:, 1)
    end

    methods
        function obj = LVPBodyTreeJAct(Joints, Bodies)
            % Call the superclass constructor
            obj = obj@BodyTree(Joints, Bodies);
            
            % Define the number of actuators
            obj.n_a = 4;

            % Define the position of the actuators in the stress-free configuration w.r.t. the base
            obj.ActuatorPosition = zeros(obj.n_a, 3);
            obj.ActuatorPosition(1, 1:3) = [0.05, 0, 0.32];
            obj.ActuatorPosition(2, 1:3) = [-0.05, 0, 0.32];
            obj.ActuatorPosition(3, 1:3) = [0, 0.05, 0.32];
            obj.ActuatorPosition(4, 1:3) = [0,-0.05, 0.32];

            % Define the body to which the actuator belongs
            obj.ActuatorBody = zeros(obj.n_a, 1);
            obj.ActuatorBody(1) = 1;
            obj.ActuatorBody(2) = 1;
            obj.ActuatorBody(3) = 1;
            obj.ActuatorBody(4) = 1;

            % Define the direction of the actuation force in the tip (local) frame
            obj.ActuatorDirection = zeros(obj.n_a, 3);
            obj.ActuatorDirection(1, 1:3) = [0, 0, -1];
            obj.ActuatorDirection(2, 1:3) = [0, 0, -1];
            obj.ActuatorDirection(3, 1:3) = [0, 0, -1];
            obj.ActuatorDirection(4, 1:3) = [0, 0, -1];

            % Find the closest centroid (index) to the actuator positions
            obj.ClosestBodyCentroid = zeros(obj.n_a, 1);
            for i = 1:obj.n_a
                bodyIdx = obj.ActuatorBody(i);
                for j = 1:obj.MaxBodiesNumber
                    if j <= obj.N_B && bodyIdx == j
                        if isnumeric(obj.Bodies{j})
                            continue;
                        end
                        Centroids = obj.Bodies{j}.Centroids;
                        [~, obj.ClosestBodyCentroid(i)] = min(vecnorm(Centroids-obj.ActuatorPosition(i, 1:3)'));
                    end
                end
            end
        end
        
        % Override the actuation matrix method
        function A = ActuationMatrix(obj, q)
            arguments (Input)
                obj (1, 1) LVPBodyTreeJAct
                q   (:, 1)
            end

            arguments (Output)
                A (:, :)
            end
                
            % Preallocate the actuation matrix
            A = zeros(obj.n, obj.n_a);

            % Update the status of the bodytree in the current configuration but update only the kinematic terms as these are the only we need
            obj.TreeUpdate(q, zeros(obj.n, 1), zeros(obj.n, 1), "EvaluateKinematicTerms", true, "EvaluateInertialTerms", false, "EvaluateExternalForces", false);

            % Iterate over all the bodies and actuators (needed for code generation)
            for i = 1:obj.n_a
                % Index of the current actuated body
                aB = obj.ActuatorBody(i);
                % Iterate over all the bodies (needed for code generation)
                for j = 1:obj.MaxBodiesNumber
                    if j <= obj.N_B && aB == j
                        if isnumeric(obj.Bodies{j})
                            continue;
                        end
                        % Get the transformation matrix
                        T = obj.Bodies{j}.T_;
                        % Get the Jacobian of the body
                        BJ = obj.Bodies{j}.BodyJacobian();
                        % Compute the effect of the actuator
                        J_Act   = T(1:3, 1:3)*BJ(1:3, :, obj.ClosestBodyCentroid(i));
                        % Compute the effect of the actuation
                        Ai      = (J_Act')*(obj.ActuatorDirection(i, 1:3)');
        
                        % Update the i-th column of the actuation matrix
                        q_start = obj.BodyConfigurationIndexes(j, 1);
                        q_end   = obj.BodyConfigurationIndexes(j, 2);
                        A(q_start:q_end, i) = Ai;
                    end
                end
            end
        end
    end
end

