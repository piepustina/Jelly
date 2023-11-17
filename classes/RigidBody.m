classdef RigidBody < Body
    %RIGIDBODY Class representing a rigid body. A rigid body is regarded as
    %a body with zero degrees of freedom. To represent a floating rigid
    %body attach a 6 degrees of freedom joint to the body.
    
    %Abstract properties implementation
    properties(Constant)
        %Number of degrees of freedom of the body.
        n = 0
    end
    properties
        % A rigid body has 10 parameters organized as follows
        % Parameters = [Mass; P_com_x; P_com_y; P_com_z; I_xx; I_yy; I_zz; I_xy; I_xz; I_yz];
        Parameters = zeros(10, 1);
    end

    methods
        function obj = RigidBody(Parameters)
            obj             = obj@Body(RigidBody.n);
            if isrow(Parameters)
                Parameters = Parameters';
            end
            obj.Parameters  = Parameters;
        end

        %Overload the methods required to describe the motion of the body
        %Center of mass position in the body frame
        function p_com_ = p_com(obj, ~)
            p_com_ = obj.Parameters(2:4);
        end
        %Inertia of the body in the body frame
        function I_ = I(obj, ~)
            I_ = [obj.Parameters(5), obj.Parameters(8) , obj.Parameters(9); 
                  obj.Parameters(8), obj.Parameters(6) , obj.Parameters(10); 
                  obj.Parameters(9), obj.Parameters(10), obj.Parameters(7)];
        end
        %Mass of the body
        function m_ = m(obj)
            m_ = obj.Parameters(1);
        end
    end
end

