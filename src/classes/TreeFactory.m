classdef TreeFactory < handle
    %"Factory" class to programmatically create joints and bodies during code generation.
    %TODO: Remove n
    %TODO: Move the methods of this class inside the TreeStructConverter.
    
    
    methods (Static)
        
        % Joint creation
        function J = CreateJoint(Type, n, Parameters)
            % Get the class constructor from the type
            Constructor = str2func(Type);
            % Call the class constructor to build the joint
            if isempty(Parameters)
                J = Constructor();
            else
                J = Constructor(Parameters);
            end
        end
        
        % Body creation
        function B = CreateBody(Type, n, Parameters)
            % Get the class contructor from the type
            Constructor = str2func(Type);
            % Call the class constructor to build the body
            if isempty(Parameters)
                B = Constructor();
            else
                B = Constructor(Parameters);
            end
        end
    end
end

