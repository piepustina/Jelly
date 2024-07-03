classdef ActuationCoordinateSystem < matlab.System
    % ActuationCoordinateSystem Implements the computation of the actuation matrix for any class extending the BodyTree and possibly implementing its own ActuationMatrix method.
    %
    
    %#codegen

    % Public, nontunable
    properties (Nontunable)
        %
        %Structure representing the robot 
        RobotStruct             = 0;
    end


    % Pre-computed constants
    properties(Access = private)
        %Internal copy of the soft robot as object
        Tree
    end

    methods
        function obj = ActuationCoordinateSystem(varargin)
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

        function y = stepImpl(obj, q)
            % Preallocate the output
            y    = cast(zeros(obj.RobotStruct.n_a, 1), 'like', q);
            
            % Compute the actuation matrix and actuator elongation
            yq   = obj.Tree.ActuationCoordinate(double(q));

            % Assign the output
            y(:) = cast(yq, 'like', q);
        end

        function resetImpl(~)
            % Initialize / reset discrete-state properties;
        end

        %% Advanced functions
        function validateInputsImpl(~, q)
            % Validate inputs to the step method at initialization
            validateattributes(q, {'single','double'},{'vector'}, 'ActuationCoordinateSystem','Configuration');
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
        
        function [out1] = getOutputSizeImpl(obj)
            %getOutputSizeImpl Return size for each output port
            out1 = [obj.RobotStruct.n_a, 1];
        end

        function [out1] = getOutputDataTypeImpl(obj)
            %getOutputDataTypeImpl Return data type for each output port
            out1 = propagatedInputDataType(obj, 1);
        end

        function [out1] = isOutputComplexImpl(obj)
            %isOutputComplexImpl Return true for each output port with complex data
            out1 = false;
        end

        function [out1] = isOutputFixedSizeImpl(obj)
            %isOutputFixedSizeImpl Return true for each output port with fixed size
            out1 = true;
        end

        function flag = isInactivePropertyImpl(~, ~)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;
        end
    end
end
