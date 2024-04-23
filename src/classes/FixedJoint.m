classdef FixedJoint < Joint
    %Class modeling a fixed joint, i.e., a hinge. 
    
    % Abstract properties implementation
    properties
        n = 0
    end
    properties
        Parameters = [];
    end
    
    methods
        % Class constructor
        function obj = FixedJoint()
            obj  = obj@Joint();
        end
    end
end

