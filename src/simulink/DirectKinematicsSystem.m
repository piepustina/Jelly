classdef DirectKinematicsSystem < matlab.System
    % DirectKinematicsSystem Implements the computation of the direct kinematics.
    %
    
    %#codegen

    % Public, nontunable
    properties (Nontunable)
        %
        % Structure representing the robot
        RobotStruct = 0;
        % Indexed of the bodies for which the direct kinematics has to be
        % evaluated
        BodyIndexes = 0;
        % Length of the body indexes
        BodyIndexesLength = 0;
    end


    % Pre-computed constants
    properties(Access = private)
        % Internal copy of the robot as object
        Tree
    end

    methods
        function obj = DirectKinematicsSystem(varargin)
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

        function T = stepImpl(obj, q)
            % Compute the direct kinematics 
            T    = cast(zeros(4*obj.BodyIndexesLength, 4), 'like', q);
            T(:) = cast(obj.Tree.DirectKinematics(double(q), obj.BodyIndexes), 'like', q);
        end

        function resetImpl(~)
            % Initialize / reset discrete-state properties;
        end

        %% Advanced functions
        function validateInputsImpl(~, q)
            % Validate inputs to the step method at initialization
            validateattributes(q, {'single','double'},{'vector'}, 'DirectKinematicsSystem','Configuration');
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
            out = [length(obj.BodyIndexes)*4, 4];
        end

        function out = getOutputDataTypeImpl(obj)
            %getOutputDataTypeImpl Return data type for each output port
            out = propagatedInputDataType(obj, 1);
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