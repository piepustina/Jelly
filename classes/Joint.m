classdef Joint < Body
    %Class that represents a generic joint of the kinematic tree.
    
    %#codegen
    methods (Access = protected)
        function obj = Joint(n)
            obj = obj@Body(n);
        end
    end

    methods (Access = public)
        %Overload the structure representation
        function s = toStruct(obj)
            % Convert object to a struct representation.
            s = struct('JointType', class(obj), 'JointParameters', obj.Parameters, 'JointDoF', obj.n);
        end
    end
end


