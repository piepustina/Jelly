classdef TreeFactory < handle
    %JOINTFACTORY Factory class that creates joints and bodies.
    %TODO: Remove n
    
    methods (Static)
        
        %Joint creation
        function J = CreateJoint(Type, n, Parameters)
            %Get the class constructor from the type
            Constructor = str2func(Type);
            %Call the class constructor to build the joint
            if isempty(Parameters{1})
                J = Constructor();
            else
                J = Constructor(Parameters{1});
            end
        end
        
        %Body creation
        function B = CreateBody(Type, n, Parameters)
            %Get the class contructor from the type
            Constructor = str2func(Type);
            %Call the class constructor to build the body
            if isempty(Parameters{1})
                B = Constructor();
            else
                B = Constructor(Parameters{1});
            end
        end
    end
end

