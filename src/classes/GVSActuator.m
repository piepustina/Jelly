classdef GVSActuator < handle
    %Abstract class modeling a thin actuator according to the Gemetric Variable Strain (GVS) theory.

    properties
        % Vector of parameters for the actuator.
        Parameters;
        % Value of the curvilinear abscissa where the actuator starts
        sStart      = 0;
        % Value of the curvilinear abscissa where the actuator ends
        sEnd        = 0;
        % Number of Gaussian points used to approximate the actuator generalized force
        NGaussPoints    = 0;
    end

    methods(Abstract)
        % Abstract method computing the distance of the actuator from the central backbone
        %
        %Args:
        %   s          (double)         : Curvilinear abscissa 
        %Return:
        %   ([3x1 - double]): Position vector describing the position of the actuator in the body local frame.
        d = ActuatorDistance(obj, s);

        % Abstract method computing the spatial derivative of the actuator distance
        %
        %Args:
        %   s          (double)         : Curvilinear abscissa 
        %Return:
        %   ([3x1 - double]): Spatial derivative of the actuator distance
        Jds = dActuatorDistance(obj, s);
        
    end
    
    methods
        function obj = GVSActuator(Parameters)
            % Construct a Geometric Variable Strain (GVS) actuator.
            %
            %Args:
            %   Parameters ([double]): Parameters of the actuator where the first two elements are :math:`s_{\mathrm{start}}, s_{\mathrm{end}}` and :math:`N_{\mathrm{Gauss}}`.
                   
            % Convert the parameters to a column vector, if needed
            if isrow(Parameters)
                Parameters = Parameters';
            end

            % Store the actuator parameters
            obj.Parameters   = Parameters;
            obj.sStart       = Parameters(1);
            obj.sEnd         = Parameters(2);
            obj.NGaussPoints = Parameters(3);
         end
        
        function s = toStruct(obj)
            %Convert the actuator object to a struct representation.
            s = struct('ActuatorType', class(obj), 'ActuatorParameters', obj.Parameters);
        end
    end

end

