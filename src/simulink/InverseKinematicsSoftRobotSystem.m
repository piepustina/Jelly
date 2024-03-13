classdef InverseKinematicsSoftRobotSystem < matlab.System
    % InverseKinematicsSoftRobotSystem Implements the computation of the inverse kinematics.
    %
    
    %#codegen

    % Public, nontunable
    properties (Nontunable)
        %
        %Structure representing the robot 
        SoftRobotStruct        = 0;
        %Points along the backbone to consider in the inverse kinematics
        BodyPoints             = 0;
        %Maximum number of iterations
        MaximumNewtonIterations = 4;
        %Task space flags
        TaskSpaceFlags          = ones(6, 1);
        %Angular error threshold
        AngularErrorThreshold   = 1e-3;
        %Linear error threshold
        LinearErrorThreshold    = 1e-2;
        %Use Gradient Descent
        UseGradientDescent (1,1) logical = false;
        %Step size for the gradient descent
        GradientDescentStepSize = 1;
        %Weight for the error
        ErrorWeight             = 1;
    end


    % Pre-computed constants
    properties(Access = private)
        %Internal copy of the robot as object
        Tree
    end

    methods
        function obj = InverseKinematicsSoftRobotSystem(varargin)
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

        function [q, converged, e] = stepImpl(obj, T, q0)
            % Preallocate the output
            q    = cast(zeros(size(q0)), 'like', q0);
            e    = zeros(6*length(obj.BodyPoints), 1);
            % Run the inverse kinematics
            [q_ik, converged, e_ik] = obj.Tree.InverseKinematics(double(T), ...
                                                                 "Points", obj.BodyPoints, ...
                                                                 "InitialGuess", double(q0), ...
                                                                 "MaxIterationNumber", obj.MaximumNewtonIterations, ...
                                                                 "TaskFlags", obj.TaskSpaceFlags, ...
                                                                 "AngularErrorThreshold", obj.AngularErrorThreshold, ...
                                                                 "LinearErrorThreshold", obj.LinearErrorThreshold, ...
                                                                 "ErrorWeight", obj.ErrorWeight, ...
                                                                 "UseGradientDescent", obj.UseGradientDescent, ...
                                                                 "GradientDescentStepSize", obj.GradientDescentStepSize);
            % Assign the output configuration
            q(:) = cast(q_ik, 'like', q0);
            e(:) = double(e_ik);
        end

        function resetImpl(~)
            % Initialize / reset discrete-state properties;
        end

        %% Advanced functions
        function flag = isInactivePropertyImpl(obj,propertyName)
            % Set properties that are not visible. Return false if property is visible based on object 
            % configuration, for the command line and System block dialog
            if strcmp(propertyName, 'GradientDescentStepSize') || strcmp(propertyName, 'ErrorWeight')
                flag = obj.UseGradientDescent == false;
            else% For all the other properties, the display flag is always ok
                flag = false;
            end
        end

        function validateInputsImpl(~, T, q0)
            % Validate inputs to the step method at initialization
            validateattributes(T , {'single','double'},    {'2d'}, 'InverseKinematicsSoftRobotSystem','Configuration');
            validateattributes(q0, {'single','double'},{'vector'}, 'InverseKinematicsSoftRobotSystem','Configuration');
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
            out1 = [obj.SoftRobotStruct.TreeStruct.n];
            out2 = 1;
            out3 = [6*length(obj.BodyPoints), 1];
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

    end
end
