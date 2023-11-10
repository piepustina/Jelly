classdef FixedJoint < JointNew
    %FIXEDJOINT Implements a fixed joint between two bodies. 
    
    %Abstract properties implementation
    properties(Constant)
        %Number of degrees of freedom of the joint
        n = 0
    end
    properties
        Parameters = [];
    end
    
    methods
        %Class constructor
        function obj = FixedJoint()
            obj  = obj@JointNew(FixedJoint.n);
        end
    end
end

