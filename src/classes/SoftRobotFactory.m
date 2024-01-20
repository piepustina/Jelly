classdef SoftRobotFactory < handle
    %"Factory" class to programmatically create actuators during code generation.
    %TODO: Move the methods of this class inside the SoftRobotStructConverter.
    
    methods (Static)
        
        % Actuator creation
        function A = CreateActuator(Type, Parameters)
            % Get the class contructor from the type
            Constructor = str2func(Type);
            % Call the class constructor to build the actuator
            if isempty(Parameters)
                A = Constructor();
            else
                A = Constructor(Parameters);
            end
        end
    end
end

