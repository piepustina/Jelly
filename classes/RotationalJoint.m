classdef RotationalJoint < JointNew
    %ROTATIONALJOINT Implements a rotational joint between two bodies. 
    
    %Abstract properties implementation
    properties(Constant)
        %Number of degrees of freedom of the joint
        n = 1
    end
    properties
        %The parameters of the joint are its DH parameters organized as
        %follows: Parameters = [alpha; a; d; theta]
        Parameters = zeros(4, 1);
    end
    
    methods
        %Class constructor
        function obj = RotationalJoint(Parameters)
            obj  = obj@JointNew(RotationalJoint.n);
            if isrow(Parameters)
                Parameters = Parameters';
            end
            obj.Parameters  = Parameters;
        end

        %Overload the methods that describe the joint motion
        %Transformation from base to body tip
        function T_ = T(obj, q)
            alpha   = obj.Parameters(1);
            a       = obj.Parameters(2);
            d       = obj.Parameters(3);
            theta   = obj.Parameters(4) + q;
            T_      =  [cos(theta), -cos(alpha)*sin(theta),  sin(alpha)*sin(theta), a*cos(theta);
                        sin(theta),  cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta);
                                 0,             sin(alpha),             cos(alpha),            d;
                                 0,                      0,                      0,            1];
        end
        %Transformation matrix from base to s
        function Ts_ = T_s(obj, q, ~)
            Ts_ = obj.T(q);
        end
        %Relative velocity of the body tip
        function v_rel_ = v_rel(obj, q, dq)
            a       = obj.Parameters(2);
            theta   = obj.Parameters(4) + q;
            dtheta  = dq;
            v_rel_ = [-a*dtheta*sin(theta);
                       a*dtheta*cos(theta);
                                         0];
        end
        %Relative angular velocity
        function omega_rel_ = omega_rel(~, ~, dq)
            omega_rel_ = [0; 0; dq];
        end
        %Relative linear acceleration of the tip 
        function a_rel_ = a_rel(obj, q, dq, ddq)
            a       = obj.Parameters(2);
            theta   = obj.Parameters(4) + q;
            dtheta  = dq;
            ddtheta = ddq;
            a_rel_  = [-a*cos(theta)*dtheta^2 - a*ddtheta*sin(theta);
                       -a*sin(theta)*dtheta^2 + a*ddtheta*cos(theta);
                                                                   0];
        end
        %Relative angular acceleration
        function domega_rel_ = domega_rel(~, ~, ~, ddq)
            domega_rel_ = [0; 0; ddq];
        end
        %Jacobian of the linear velocity of the tip with respect to q in the tip frame
        function v_par_ = v_par(obj, ~)
            alpha   = obj.Parameters(1);
            a       = obj.Parameters(2);
            v_par_ = [           0;
                      a*cos(alpha);
                     -a*sin(alpha)];
        end
        %Jacobian of the angular velocity with respect to q in the tip frame
        function omega_par_ = omega_par(obj, ~)
            alpha       = obj.Parameters(1);
            omega_par_  = [0; sin(alpha); cos(alpha)];
        end
    end
end

