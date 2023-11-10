classdef TreeFactory < handle
    %JOINTFACTORY Factory class that creates joints and bodies of the various type
    
    methods (Static)
        
        function J = CreateJoint(Type, n, Parameters)
            %Get the class constructor from the type
            Constructor = str2func(Type);
            %Call the class constructor
            J = Constructor(Parameters{1});
            % switch Type
            %     case prova.PCCCylindricalDelta_Joint.Type
            %         J = prova.PCCCylindricalDelta_Joint(Parameters{1}(1));
            %     case prova.PCCPlanar_Joint.Type
            %         J = prova.PCCPlanar_Joint(Parameters{1}(1));
            %     case prova.PCCPlanar_Joint_no_elongation.Type
            %         J = prova.PCCPlanar_Joint_no_elongation(Parameters{1}(1));
            %     case prova.PACCylindrical_Delta_Joint.Type
            %         J = prova.PACCylindrical_Delta_Joint(Parameters{1}(1));
            %     case prova.Rotational_Joint.Type
            %         J = prova.Rotational_Joint(Parameters{1});
            %     case prova.PCSCylindrical_Joint.Type
            %         J = prova.PCSCylindrical_Joint(Parameters{1}(1));
            %     case prova.GVSJoint.Type
            %         J = prova.GVSJoint(Parameters{2}, n, Parameters{1});
            %     otherwise 
            %         error("Joint not yet implemented");
            % end
        end

        function B = CreateBody(Type, n, Parameters)
            Constructor = str2func(Type);
            B = Constructor(Parameters{1});
            % switch Type
            %     case prova.PCCCylindricalDelta_Body.Type
            %         B = prova.PCCCylindricalDelta_Body(Parameters{1});
            %     case prova.PCCPlanar_Body.Type
            %         B = prova.PCCPlanar_Body(Parameters{1});
            %     case prova.PCCPlanar_Body_no_elongation.Type
            %         B = prova.PCCPlanar_Body_no_elongation(Parameters{1});
            %     case prova.PACCylindrical_Delta_Body.Type
            %         B = prova.PACCylindrical_Delta_Body(Parameters{1});
            %     case prova.Rotational_Body.Type
            %         B = prova.Rotational_Body(Parameters{1});
            %     case prova.PCSCylindrical_Body.Type
            %         B = prova.PCSCylindrical_Body(Parameters{1});
            %     case prova.GVSBody.Type
            %         B = prova.GVSBody(Parameters{2}, n, Parameters{1});
            %     otherwise 
            %         error("Body not yet implemented");
            % end
        end
    end
end

