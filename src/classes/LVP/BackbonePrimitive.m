classdef (Abstract) BackbonePrimitive < handle 
    % Abstract class that is used to model if a primitive acts on the backbone or not
    
    methods
        % Strain basis and its derivative w.r.t. s
        [J, dJ]  = StrainBasis(obj, s);
    end

end

