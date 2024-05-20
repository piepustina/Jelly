classdef LVPBodyTree < BodyTree
    %LVPBODYTREE Test class used for defining custom actuation matrix
    
    properties (Access = private)
        ActuatorPosition  (:, 3)
        ActuatorBody      (:, 1)
        ActuatorDirection (:, 3)
    end

    methods
        function obj = LVPBodyTree(Joints, Bodies)
            % Call the superclass constructor
            obj = obj@BodyTree(Joints, Bodies);
            
            % Define the number of actuators
            obj.n_a = 4;

            % Define the position of the actuators in the stress-free configuration and tip (local) frame
            obj.ActuatorPosition = zeros(obj.n_a, 3);
            obj.ActuatorPosition(1, 1:3) = [0.015, 0, 0];
            obj.ActuatorPosition(2, 1:3) = [-0.015, 0, 0];
            obj.ActuatorPosition(3, 1:3) = [0, 0.015, 0];
            obj.ActuatorPosition(4, 1:3) = [0,-0.015, 0];

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
        end
        
        % Override the actuation matrix method
        function A = ActuationMatrix(obj, q)
            arguments (Input)
                obj (1, 1) LVPBodyTree
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
                        % Linear Jacobian of the body tip in the body frame
                        J_L     = obj.Bodies{j}.v_par_;
                        % Angular Jacobian of the body tip in the body frame
                        J_A     = obj.Bodies{j}.omega_par_;
                        % Linear Jacobian for the actuator
                        J_Act   = J_L - skew([1; 0; 0])*J_A*obj.ActuatorPosition(i, 1)...
                                      - skew([0; 1; 0])*J_A*obj.ActuatorPosition(i, 2)...
                                      - skew([0; 0; 1])*J_A*obj.ActuatorPosition(i, 3);
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

