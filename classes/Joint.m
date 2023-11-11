classdef Joint < Body
    %JOINT Class that represents a generic Joint. The class is represented
    %as a body with no inertial parameters.
    %#codegen

    methods (Access = protected)
        function obj = Joint(n)
            obj = obj@Body(n);
        end
    end

    methods (Access = public)
        %Overload the structure representation
        function s = toStruct(obj)
            s = struct('JointType', class(obj), 'JointParameters', {{obj.Parameters}}, 'JointDoF', obj.n);
        end
    end
end


