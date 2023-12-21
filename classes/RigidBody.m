classdef RigidBody < Body
    %Class modeling a rigid body. 

    % A rigid body can be modeled as a Body class object with zero DoF. 
    % To represent a floating rigid body, attach a 6 degrees of freedom joint to the body.
    
    % Abstract properties implementation
    properties(Constant)
        n = 0
    end
    properties
        %Vector collecting the parameters of the rigid body as :math:`\mathrm{Parameters} = \left( m \,\, p_{\mathrm{CoM}_{x}} \,\, p_{\mathrm{CoM}_{y}} \,\, p_{\mathrm{CoM}_{z}} \,\, I_{xx} \,\, I_{yy} \,\, I_{zz} \,\, I_{xy} \,\, I_{xz} \,\, I_{yz} \right)^{T} \in \mathbb{R}^{10 \times 1}`
        %where :math:`m` is the body mass, :math:`p_{\mathrm{CoM}_{x}}, p_{\mathrm{CoM}_{y}}` and :math:`p_{\mathrm{CoM}_{z}}` denote the center of mass position in the body distal end frame, and :math:`I_{xx}, I_{yy},
        %I_{zz}, I_{xy}, I_{xz}` and :math:`I_{yz}` are the components of the body inertia matrix expressed in the body distal end frame.
        Parameters = zeros(10, 1);
    end

    methods
        function obj = RigidBody(Parameters)
            %Construct a rigid body.
            %
            %Args:
            %   Parameters ([double], [sym]): Parameters of the body, specified as :math:`m, p_{\mathrm{CoM}_{x}}, p_{\mathrm{CoM}_{y}}, p_{\mathrm{CoM}_{z}}, I_{xx}, I_{yy}, I_{zz}, I_{xy}, I_{xz}` and :math:`I_{yz}`

            obj             = obj@Body(RigidBody.n);
            if isrow(Parameters)
                Parameters = Parameters';
            end
            obj.Parameters  = Parameters;
        end

        %% Overload the methods required to describe the motion of the body
        % Center of mass position in the body frame
        function p_com_ = p_com(obj, ~)
            p_com_ = obj.Parameters(2:4);
        end
        % Inertia of the body in the body frame
        function I_ = I(obj, ~)
            I_ = [obj.Parameters(5), obj.Parameters(8) , obj.Parameters(9); 
                  obj.Parameters(8), obj.Parameters(6) , obj.Parameters(10); 
                  obj.Parameters(9), obj.Parameters(10), obj.Parameters(7)];
        end
        % Mass of the body
        function m_ = m(obj)
            m_ = obj.Parameters(1);
        end
    end
end

