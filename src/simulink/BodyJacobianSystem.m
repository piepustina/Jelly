classdef BodyJacobianSystem < matlab.System
    % BodyJacobianSystem Implements the computation of the elastic force.
    %
    
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
        function obj = BodyJacobianSystem(varargin)
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

        function J = stepImpl(obj, q)
            % Compute the body Jacobian
            J    = cast(zeros(6*obj.Tree.N_B, length(q)), 'like', q);
            J(:) = cast(obj.Tree.BodyJacobian(double(q)), 'like', q);
        end

        function resetImpl(~)
            % Initialize / reset discrete-state properties;
        end

        %% Advanced functions
        function validateInputsImpl(~, q)
            % Validate inputs to the step method at initialization
            validateattributes(q, {'single','double'},{'vector'}, 'BodyJacobianSystem','Configuration');
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
            out = [obj.RobotStruct.n];
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
