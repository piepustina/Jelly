classdef MassInverseDynamicsSystem < matlab.System
    % InverseDynamics Implements the inverse dynamics algorithm.
    %#codegen

    % Public, nontunable
    properties (Nontunable)
        %
        %robot_struct struct = struct('Mass', 1, 'Type', 'Rigid');
        RobotStruct = 0;
        % Evaluate generalized external forces such as damping and elastic forces.
        EvaluateExternalForces (1, 1) logical = true
    end


    % Pre-computed constants
    properties(Access = private)
        %Internal copy of the robot as object
        Tree
    end

    methods
        function obj = MassInverseDynamicsSystem(varargin)
            setProperties(obj, nargin, varargin{:})
        end
    end


    methods(Access = protected)
      

        function setupImpl(obj)
            obj.Tree = TreeStructConverter.StructToObject(obj.RobotStruct);
        end

        function [M, tau] = stepImpl(obj, q, dq, ddq)
            % Output preallocation
            M           = cast(zeros(numel(q), numel(q)), 'like', q);
            tau         = cast(zeros(size(q)), 'like', q);
            
            % Run a step of the inverse dynamics and compute both the mass matrix and the generalized actuation forces 
            [Mq, tauq]  = obj.Tree.UnifiedInverseDynamics(double(q), double(dq), double(ddq), "EvaluateExternalForces", obj.EvaluateExternalForces);
            M(:, :)     = cast(Mq, 'like', q);
            tau(:)      = cast(tauq, 'like', q);
            
        end

        function resetImpl(~)
            % Initialize / reset discrete-state properties;
        end

        %% Advanced functions
        function validateInputsImpl(~, q, dq, ddq)
            % Validate inputs to the step method at initialization
            validateattributes(q, {'single','double'},{'vector'},'MassInverseDynamicsSystem','Configuration');
            validateattributes(dq,{'single','double'},{'vector'},'MassInverseDynamicsSystem','Configuration Velocity');
            validateattributes(ddq, {'single','double'},{'vector'},'MassInverseDynamicsSystem','Configuration Acceleration');
        end

        function validatePropertiesImpl(~)
            % Validate related or interdependent property values
        end

        function ds = getDiscreteStateImpl(~)
            % Return structure of properties with DiscreteState attribute
            ds = struct([]);
        end

        function processTunedPropertiesImpl(~)
            % Perform actions when tunable properties change
            % between calls to the System object
        end

        function flag = isInputSizeMutableImpl(~,~)
            %isInputSizeMutableImpl Specify input size mutability
            %   Return false if input size cannot change
            %   between calls to the System object
            flag = false;
        end

        function flag = isInputDataTypeMutableImpl(~,~)
            %isInputDataTypeMutableImpl Specify input type mutability
            %   Return false if input data type cannot change
            %   between calls to the System object
            flag = false;
        end

        function num = getNumOutputsImpl(obj)
            %getNumOutputsImpl Define total number of outputs
            num = 2;
        end
        
        function [out1, out2] = getOutputSizeImpl(obj)
            %getOutputSizeImpl Return size for each output port
            out1 = [obj.RobotStruct.n obj.RobotStruct.n];
            out2 = [obj.RobotStruct.n 1];
        end

        function [out1, out2] = getOutputDataTypeImpl(obj)
            %getOutputDataTypeImpl Return data type for each output port
            out1 = propagatedInputDataType(obj, 1);
            out2 = propagatedInputDataType(obj, 1);

        end

        function [out1, out2] = isOutputComplexImpl(obj)
            %isOutputComplexImpl Return true for each output port with complex data
            out1 = false;
            out2 = false;
            
        end

        function [out1, out2] = isOutputFixedSizeImpl(obj)
            %isOutputFixedSizeImpl Return true for each output port with fixed size
            out1 = true;
            out2 = true;
        end

        function flag = isInactivePropertyImpl(~, ~)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;
        end
    end
end
