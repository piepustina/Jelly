classdef InverseKinematicsSystem < matlab.System
    % InverseKinematicsSystem Implements the computation of the inverse kinematics.
    %
    
    %#codegen

    % Public, nontunable
    properties (Nontunable)
        %
        %Structure representing the robot 
        RobotStruct             = 0;
        %Indexes of the bodies to consider in the inverse kinematics
        BodyIndexes             = 1;
        %Maximum number of Newton iterations
        MaximumNewtonIterations = 4;
        %Task space flags
        TaskSpaceFlags          = ones(6, 1);
    end


    % Pre-computed constants
    properties(Access = private)
        %Internal copy of the robot as object
        Tree
    end

    methods
        function obj = InverseKinematicsSystem(varargin)
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

        function [q, converged, e] = stepImpl(obj, T, q0)
            % Preallocate the output
            q    = cast(zeros(size(q0)), 'like', q0);
            e    = zeros(6*length(obj.BodyIndexes), 1);
            % Run the inverse kinematics
            [q_ik, converged, e_ik] = obj.Tree.InverseKinematics(double(T), ...
                                                "BodyIndexes", obj.BodyIndexes, ...
                                                "InitialGuess", double(q0), ...
                                                "MaxIterationNumber", obj.MaximumNewtonIterations, ...
                                                "TaskFlags", obj.TaskSpaceFlags);
            % Assign the output configuration
            q(:) = cast(q_ik, 'like', q0);
            e(:) = double(e_ik);
        end

        function resetImpl(~)
            % Initialize / reset discrete-state properties;
        end

        %% Advanced functions
        function validateInputsImpl(~, T, q0)
            % Validate inputs to the step method at initialization
            validateattributes(T , {'single','double'},    {'2d'}, 'InverseKinematicsSystem','Configuration');
            validateattributes(q0, {'single','double'},{'vector'}, 'InverseKinematicsSystem','Configuration');
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
            num = 3;
        end
        
        function [out1, out2, out3] = getOutputSizeImpl(obj)
            %getOutputSizeImpl Return size for each output port
            out1 = [obj.RobotStruct.n];
            out2 = 1;
            out3 = [6*length(obj.BodyIndexes), 1];
        end

        function [out1, out2, out3] = getOutputDataTypeImpl(obj)
            %getOutputDataTypeImpl Return data type for each output port
            out1 = propagatedInputDataType(obj, 2);
            out2 = "double";
            out3 = "double";
        end

        function [out1, out2, out3] = isOutputComplexImpl(obj)
            %isOutputComplexImpl Return true for each output port with complex data
            out1 = false;
            out2 = false;
            out3 = false;
        end

        function [out1, out2, out3] = isOutputFixedSizeImpl(obj)
            %isOutputFixedSizeImpl Return true for each output port with fixed size
            out1 = true;
            out2 = true;
            out3 = true;
        end

        function flag = isInactivePropertyImpl(~, ~)
            % Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            flag = false;
        end
    end
end
