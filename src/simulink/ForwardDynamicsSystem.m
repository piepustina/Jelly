classdef ForwardDynamicsSystem < matlab.System
    % Forward Dynamics Implements the forward dynamics algorithm for a soft
    % robot.
    %
    % This template includes the minimum set of functions required
    % to define a System object with discrete state.
    
    %#codegen

    % Public, nontunable
    properties (Nontunable)
        %
        %robot_struct struct = struct('Mass', 1, 'Type', 'Rigid');
        RobotStruct = 0;
    end


    % Pre-computed constants
    properties(Access = private)
        %Internal copy of the robot as object
        Tree
    end

    methods
        function obj = ForwardDynamicsSystem(varargin)
            %ForwardDynamics Constructor for Forward Dynamics block system object
            
            % Support name-value pair arguments when constructing object
            setProperties(obj, nargin, varargin{:})
        end
    end


    methods(Access = protected)
      

        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            %
            %Convert the robot into a Tree
            obj.Tree = TreeStructConverter.StructToObject(obj.RobotStruct);
        end

        function ddq = stepImpl(obj, q, dq, u)
            % Run a step of the forward dynamics
            ddq = cast(zeros(size(q)), 'like', q);
            ddq(:) = cast(obj.Tree.ForwardDynamics(double(q), double(dq), double(u)), 'like', q);
        end

        function resetImpl(~)
            % Initialize / reset discrete-state properties;
        end

        %% Advanced functions
        function validateInputsImpl(~, q, dq, u)
            % Validate inputs to the step method at initialization
            validateattributes(q, {'single','double'},{'vector'},'ForwardDynamicsSystem','Configuration');
            validateattributes(dq,{'single','double'},{'vector'},'ForwardDynamicsSystem','Configuration Velocity');
            validateattributes(u, {'single','double'},{'vector'},'ForwardDynamicsSystem','Input Force');
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
            num = 1;
        end
        
        function out = getOutputSizeImpl(obj)
            %getOutputSizeImpl Return size for each output port
            out = [obj.RobotStruct.n 1];
        end

        function out = getOutputDataTypeImpl(obj)
            %getOutputDataTypeImpl Return data type for each output port
            out = propagatedInputDataType(obj,1);
        end

        function out = isOutputComplexImpl(obj)
            %isOutputComplexImpl Return true for each output port with complex data
            out = false;
        end

        function out = isOutputFixedSizeImpl(obj)
            %isOutputFixedSizeImpl Return true for each output port with fixed size
            out = true;
        end

        function flag = isInactivePropertyImpl(~, ~)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;
        end
    end
end
