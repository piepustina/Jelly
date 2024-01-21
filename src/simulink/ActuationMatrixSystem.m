classdef ActuationMatrixSystem < matlab.System
    % ActuationMatrixSystem Implements the computation of the actuation matrix for a soft robot.
    %
    
    %#codegen

    % Public, nontunable
    properties (Nontunable)
        %
        %Structure representing the robot 
        SoftRobotStruct             = 0;
    end


    % Pre-computed constants
    properties(Access = private)
        %Internal copy of the soft robot as object
        Tree
    end

    methods
        function obj = ActuationMatrixSystem(varargin)
            setProperties(obj, nargin, varargin{:})
        end
    end


    methods(Access = protected)
      

        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            %
            %Convert the robot into a Tree
            obj.Tree = SoftRobotStructConverter.StructToObject(obj.SoftRobotStruct);
        end

        function [A, y] = stepImpl(obj, q)
            % Preallocate the output
            A    = cast(zeros(obj.SoftRobotStruct.TreeStruct.n, obj.SoftRobotStruct.N_A), 'like', q);
            y    = cast(zeros(obj.SoftRobotStruct.N_A, 1), 'like', q);
            % Compute the actuation matrix and actuator elongation
            [Aq, yq] = obj.Tree.ActuationMatrix(double(q));

            %disp(Aq);
            %disp(yq);

            % Assign the output
            A(:) = cast(Aq, 'like', q);
            y(:) = cast(yq, 'like', q);
        end

        function resetImpl(~)
            % Initialize / reset discrete-state properties;
        end

        %% Advanced functions
        function validateInputsImpl(~, q)
            % Validate inputs to the step method at initialization
            validateattributes(q, {'single','double'},{'vector'}, 'ActuationMatrixSystem','Configuration');
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
            out1 = [obj.SoftRobotStruct.TreeStruct.n, obj.SoftRobotStruct.N_A];
            out2 = [obj.SoftRobotStruct.N_A];
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
