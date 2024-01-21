classdef SoftRobotStructConverter < handle
    %Utility class to convert a SoftRobot object into a struct and viceversa. 
    %This class is used for code generation in Simulink.
    
    %#codegen
    
    % Define the maximum size of the parameters of the struct (required for code generation). 
    properties (Constant)
        StringSize      = 40;% Used for the ActuatorType and JointType fields
        ParametersSize  = 20;% Used for the size of the parameters
    end
    
    methods (Static)
        function SoftRobotStruct = ObjectToStruct(R)
            % R is an instance of the class SoftRobot or its a class inheriting from it
            if ~isa(R, 'SoftRobot')
                error("R is not an instance of SoftRobot.m");
            end
            % Get the struct representation of the robot as a BodyTree
            BodyTreeStruct = TreeStructConverter.ObjectToStruct(R);

            % Add also a representation of the actuators
            ActuatorStruct = struct('ActuatorType', char('#'*ones(1, SoftRobotStructConverter.StringSize)),...
                                    'LengthActuatorType', 0, ...
                                    'ActuatorParameters', nan(SoftRobotStructConverter.ParametersSize, 1), ...
                                    'NumberActuatorParameters', 0);
            % Structure representation of the soft robot
            SoftRobotStruct = struct('TreeStruct', BodyTreeStruct, ...
                                     'Actuators', repmat(ActuatorStruct, R.N_A, 1), ...
                                     'N_A', R.N_A);
            % Iterate over the bodies of the tree
            for i = 1:R.N_A
                % Get the actuator representation as a struct
                S_actuator  = R.Actuators{i}.toStruct();
                
                % Save the actuator data
                cActuatorType = char(S_actuator.ActuatorType);
                SoftRobotStruct.Actuators(i).LengthActuatorType = length(cActuatorType);
                SoftRobotStruct.Actuators(i).ActuatorType(1:SoftRobotStruct.Actuators(i).LengthActuatorType) = cActuatorType;
                SoftRobotStruct.Actuators(i).NumberActuatorParameters = length(S_actuator.ActuatorParameters);
                SoftRobotStruct.Actuators(i).ActuatorParameters(1:SoftRobotStruct.Actuators(i).NumberActuatorParameters) = S_actuator.ActuatorParameters;
            end
        end
        
        % Converts a structure representing a SoftRobot into the corresponding
        % object. To comply with the code generation it creates additional fake
        % Actuators to fill the maximum number allowed. 
        function SoftRobotObject = StructToObject(S)
            % Allocate the actuators
            Actuators = cell(SoftRobot.MaxActuatorsNumber, 1);
            
            N_A = S.N_A;
            for i = 1:N_A
                % Create the actuator
                ActuatorType = S.Actuators(i).ActuatorType(1:S.Actuators(i).LengthActuatorType);
                ActuatorParameters = S.Actuators(i).ActuatorParameters(1:S.Actuators(i).NumberActuatorParameters);
                Actuators{i} = SoftRobotFactory.CreateActuator(ActuatorType, ActuatorParameters);
            end
            
            % Fill the remainig actuators with zeros
            for i = N_A+1:SoftRobot.MaxActuatorsNumber
                Actuators{i} = 0;
            end

            % Construct the bodytree object
            BodyTreeObject                      = TreeStructConverter.StructToObject(S.TreeStruct);
            
            % Construct the soft robot object by adding also the actuators
            SoftRobotObject                     = SoftRobot(BodyTreeObject.Joints, BodyTreeObject.Bodies, Actuators);
            SoftRobotObject.g                   = BodyTreeObject.g;
            SoftRobotObject.T0                  = BodyTreeObject.T0;
        end
    end
end

