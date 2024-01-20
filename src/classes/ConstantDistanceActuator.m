classdef ConstantDistanceActuator < GVSActuator
    %Class modeling a thin actuator at a constant distance from the robot backbone.
    
    methods
        function obj = ConstantDistanceActuator(Parameters)
            %Construct a actuator with a constant distance from the robot backbone.
            %
            %Args:
            %   Parameters ([double]): Parameters of the actuator, specified as :math:`s_{\mathrm{start}}, s_{\mathrm{end}}, N_{\mathrm{Gauss}}` and :math:`d`.
            
            obj = obj@GVSActuator(Parameters);
        end
        
        % Implement abstract methods
        function d = ActuatorDistance(obj, s)
            sStart = obj.Parameters(1); sEnd = obj.Parameters(2);
            if s >= sStart && s <= sEnd
                d = obj.Parameters(4:6);
            else
                d = zeros(3, 1);
            end
        end

        function dd = dActuatorDistance(obj, s)
            dd = zeros(3, 1);
        end
    end
end

