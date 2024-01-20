classdef TreeStructConverter < handle
    %Utility class to convert a BodyTree object into a struct and viceversa. 
    %This class is used for code generation in Simulink.
    
    %#codegen
    
    % Define the maximum size of the parameters of the struct (required for code generation). 
    properties (Constant)
        StringSize      = 40;% Used for the BodyType and JointType fields
        ParametersSize  = 20;% Used for the size of the parameters
    end
    
    methods (Static)
        function TreeStruct = ObjectToStruct(R)
            % R is an instance of the class BodyTree or its a class inheriting from it
            if ~isa(R, 'BodyTree')
                error("R is not an instance of BodyTree.m");
            end
            % Iterate over the bodies of the robot
            BodyStruct = struct('BodyType', char('#'*ones(1, TreeStructConverter.StringSize)),...
                                'LengthBodyType', 0, ...
                                'BodyParameters', nan(TreeStructConverter.ParametersSize, 1), ...
                                'NumberBodyParameters', 0, ...
                                'BodyDoF', 0, ...
                                'JointType', char('#'*ones(1, TreeStructConverter.StringSize)), ...
                                'LengthJointType', 0, ...
                                'JointParameters', nan(TreeStructConverter.ParametersSize, 1),...
                                'NumberJointParameters', 0, ...
                                'JointDoF', 0);
            TreeStruct = struct('n', R.n, ...
                                'g', R.g, ...
                                'T0', R.T0, ...
                                'MassConditionNumber', R.MassConditionNumber, ...
                                'Bodies', repmat(BodyStruct, R.N_B, 1));
            % Iterate over the bodies of the tree
            for i = 1:R.N_B
                % Get the body and joint representation as structs
                S_body  = R.Bodies{i}.toStruct();
                S_joint = R.Joints{i}.toStruct();
                
                % Save the body
                cBodyType = char(S_body.BodyType);
                TreeStruct.Bodies(i).LengthBodyType = length(cBodyType);
                TreeStruct.Bodies(i).BodyType(1:TreeStruct.Bodies(i).LengthBodyType) = cBodyType;
                TreeStruct.Bodies(i).NumberBodyParameters = length(S_body.BodyParameters);
                TreeStruct.Bodies(i).BodyParameters(1:TreeStruct.Bodies(i).NumberBodyParameters) = S_body.BodyParameters;
                TreeStruct.Bodies(i).BodyDoF = S_body.BodyDoF;
                % Save the joint
                cJointType = char(S_joint.JointType);
                TreeStruct.Bodies(i).LengthJointType = length(cJointType);
                TreeStruct.Bodies(i).JointType(1:TreeStruct.Bodies(i).LengthJointType) = cJointType;
                TreeStruct.Bodies(i).NumberJointParameters = length(S_joint.JointParameters);
                TreeStruct.Bodies(i).JointParameters(1:TreeStruct.Bodies(i).NumberJointParameters) = S_joint.JointParameters;
                TreeStruct.Bodies(i).JointDoF = S_joint.JointDoF;
            end
        end
        
        % Converts a structure representing a Tree into the corresponding
        % object. To comply with the code generation it creates additional
        % Bodies and Joints to fill the maximum number allowed. 
        function TreeObject = StructToObject(S)
            Joints = cell(BodyTree.MaxBodiesNumber, 1);
            Bodies = cell(BodyTree.MaxBodiesNumber, 1);
            
            N_B = length(S.Bodies);
            for i = 1:N_B
                % Create the joint
                JointType = S.Bodies(i).JointType(1:S.Bodies(i).LengthJointType);
                JointParameters = S.Bodies(i).JointParameters(1:S.Bodies(i).NumberJointParameters);
                Joints{i} = TreeFactory.CreateJoint(JointType, S.Bodies(i).JointDoF, JointParameters);
                % Create the body
                BodyType = S.Bodies(i).BodyType(1:S.Bodies(i).LengthBodyType);
                BodyParameters = S.Bodies(i).BodyParameters(1:S.Bodies(i).NumberBodyParameters);
                Bodies{i} = TreeFactory.CreateBody(BodyType, S.Bodies(i).BodyDoF, BodyParameters);
            end
            
            % Fill the remainig bodies with zeros
            for i = N_B+1:BodyTree.MaxBodiesNumber
                Bodies{i} = 0;
                Joints{i} = 0;
            end

            TreeObject = BodyTree(Joints, Bodies);
            TreeObject.g = S.g;
            TreeObject.T0 = S.T0;
            TreeObject.MassConditionNumber = S.MassConditionNumber;
        end
    end
end

