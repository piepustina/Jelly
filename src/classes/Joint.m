classdef Joint < Body
    %Abstract class modeling a generic joint of a :class:`BodyTree`. 
    
    %#codegen
    methods (Access = protected)
        function obj = Joint()
            %Construct the joint.
            %
            %Args:
            %   n (double): Number of DoF of the joint
            obj = obj@Body();
        end
    end

    methods (Access = public)
        % Overload the structure representation
        function s = toStruct(obj)
            %Convert the joint object to a struct representation.
            s = struct('JointType', class(obj), 'JointParameters', obj.Parameters, 'JointDoF', obj.n);
        end
    end
end


