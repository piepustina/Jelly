classdef ConstantDistanceActuator < GVSActuator
    %Class modeling a thin actuator at a constant distance from the robot backbone.
    
    methods
        function obj = ConstantDistanceActuator(Parameters)
            %Construct a actuator with a constant distance from the robot backbone.
            %
            %Args:
            %   Parameters ([double]): Parameters of the actuator, specified as :math:`L_0, N_{\mathrm{Gauss}}` and :math:`d`.
            
            obj = obj@GVSActuator(Parameters);
        end
        
        % Implement abstract methods
        function d = ActuatorDistance(obj, s)
            d = obj.Parameters(3:5);
        end

        function dd = dActuatorDistance(obj, s)
            dd = zeros(3, 1);
        end
    end
end

