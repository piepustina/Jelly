classdef LVPBodyTree < BodyTree
    %LVPBODYTREE Test class used for defining custom actuation matrix
    
    methods
        function obj = LVPBodyTree(Joints, Bodies)
            % Call the superclass constructor
            obj = obj@BodyTree(Joints, Bodies);
        end
        
        % Override the actuation matrix method
        function A = ActuationMatrix(obj, q)
            arguments (Input)
                obj (1, 1) LVPBodyTree
                q   (:, 1)
            end

            arguments (Output)
                A (:, :)
            end
                
            disp("Metodo sovrascritto!");
            A = eye(obj.n);
        end
    end
end

