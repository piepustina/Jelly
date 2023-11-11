classdef Body < handle
    %BODY 
    properties (Access = public)
        %Kinetic terms
        T_                  = eye(4)
        v_rel_              = zeros(3, 1)
        omega_rel_          = zeros(3, 1)
        a_rel_              = zeros(3, 1)
        domega_rel_         = zeros(3, 1)
        v_par_
        omega_par_
        %Inertial terms
        p_com_              = zeros(3, 1)
        v_com_rel_          = zeros(3, 1)
        a_com_rel_          = zeros(3, 1)
        m_                  = 0
        I_                  = zeros(3, 3)
        J_                  = zeros(3, 3)
        int_dr_             = zeros(3, 1)
        int_ddr_            = zeros(3, 1)
        int_r_X_dr_         = zeros(3, 1)
        int_r_X_ddr_        = zeros(3, 1)
        int_dr_X_pv_r_
        int_pv_r_O_dd_r_
        int_dr_O_dr_        = 0
        grad_int_dr_
        grad_int_r_X_dr_
        grad_J_
        grad_v_com_
        xi_                 = zeros(6, 1);
        K_
        D_
    end

    properties(Abstract, Constant)
        %Number of degrees of freedom of the body
        n;
    end

    properties(Abstract)
        %Parameters of the body
        Parameters;
    end

    %Methods required to update the body state
    methods(Access = public)
            %% Kinematic terms
            %Transformation matrix from tip to the base
            function T_ = T(~, q)
                T_ = eye(4, 'like', q);
            end
            %Transformation matrix from s to the base
            function T_s_ = T_s(~, q, ~)
                T_s_ = eye(4, 'like', q);
            end
            %Relative velocity of the tip in the base frame
            function v_rel_ = v_rel(~, q, ~)
                v_rel_ = zeros(3, 1, 'like', q);
            end
            %Relative angular velocity of the body in the base frame
            function omega_rel_ = omega_rel(~, q, ~)
                omega_rel_ = zeros(3, 1, 'like', q);
            end
            %Relative acceleration of the tip in the base frame
            function a_rel_ = a_rel(~, q, ~, ~)
                a_rel_ = zeros(3, 1, 'like', q);
            end
            %Relative angular acceleration of the body
            function domega_rel_ = domega_rel(~, q, ~, ~)
                domega_rel_ = zeros(3, 1, 'like', q);
            end
            %Jacobian of the linear tip velocity in the tip frame
            function v_par_ = v_par(obj, q)
                v_par_ = zeros(3, obj.n, 'like', q);
            end
            %Jacobian of the angular velocity in the tip frame
            function omega_par_ = omega_par(obj, q)
                omega_par_ = zeros(3, obj.n, 'like', q);
            end
            
            %% Inertial terms
            %Center of mass position in the body frame
            function p_com_ = p_com(~, q)
                p_com_ = zeros(3, 1, 'like', q);
            end
            %Linear velocity of the center of mass in the body frame
            function v_com_rel_ = v_com_rel(~, q, ~)
                v_com_rel_ = zeros(3, 1, 'like', q);
            end
            %Linear acceleration of the center of mass in the body frame
            function a_com_rel_ = a_com_rel(~, q, ~, ~)
                a_com_rel_ = zeros(3, 1, 'like', q);
            end
            %Inertia of the body in the body frame
            function I_ = I(~, q)
                I_ = zeros(3, 3, 'like', q);
            end
            %Mass of the body
            function m_ = m(~)
                m_ = 0;
            end
            %Time derivative of the inertia in the body frame
            function J_ = J(~, q, ~)
                J_ = zeros(3, 3, 'like', q);
            end
            %Integral of \dot{r}
            function int_dr_ = int_dr(~, q, ~)
                %TODO: This is always zero, remove
                int_dr_ = zeros(3, 1, 'like', q);
            end
            %Integral of \ddot{r}
            function int_ddr_ = int_ddr(~, q, ~, ~)
                %TODO: This is always zero, remove
                int_ddr_ = zeros(3, 1, 'like', q);
            end
            %Integral of \cross(r, \dot{r})
            function int_r_X_dr_ = int_r_X_dr(~, q, ~)
                int_r_X_dr_ = zeros(3, 1, 'like', q);
            end
            %Integral of \cross(r, \ddot{r})
            function int_r_X_ddr_ = int_r_X_ddr(~, q, ~, ~)
                int_r_X_ddr_ = zeros(3, 1, 'like', q);
            end
            %Integral of \cross(\dor{r}, \jacobian{r}{q})
            function int_dr_X_pv_r_ = int_dr_X_pv_r(obj, q, ~)
                int_dr_X_pv_r_ = zeros(obj.n, 3, 'like', q);
            end
            %Integral of \dot(\jacobian{r}{q}, \ddot{r})
            function int_pv_r_O_dd_r_ = int_pv_r_O_dd_r(obj, q, ~, ~)
                int_pv_r_O_dd_r_ = zeros(obj.n, 1, 'like', q);
            end
            %Integral of \dot(\dot{r}, \dot{r})
            function int_dr_O_dr_ = int_dr_O_dr(~, ~, ~)
                %TODO: This is always zero, remove
                int_dr_O_dr_ = 0;
            end
            %Jacobian of the integral of \dot{r}
            function grad_int_dr_ = grad_int_dr(obj, q)
                %TODO: This is always zero, remove
                grad_int_dr_ = zeros(obj.n, 3, 'like', q);
            end
            %Jacobian of the integral of \cross{r, \dot{r}}
            function grad_int_r_X_dr_ = grad_int_r_X_dr(obj, q)
                grad_int_r_X_dr_ = zeros(obj.n, 3, 'like', q);
            end
            %Jacobian of the time derivative of the inertia
            function grad_J_ = grad_J(obj, q)
                grad_J_ = zeros(3, 3, obj.n, 'like', q);
            end
            %Jacobian of the center of mass velocity
            function grad_v_com_ = grad_v_com(obj, q)
                grad_v_com_ = zeros(obj.n, 3, 'like', q);
            end
            %Strain function
            function xi_ = xi(~, q, ~)
                xi_ = zeros(6, 1, 'like', q);
            end
            %Generalized stiffness force
            function K_ = K(obj, q)
                K_ = zeros(obj.n, 1, 'like', q);
            end
            %Generalized damping force
            function D_ = D(obj, q, ~)
                D_ = zeros(obj.n, 1, 'like', q);
            end
    end

    methods (Access = protected)
        %Class constructor
        function obj = Body(n)
            %Initialize all the variables
            obj.v_par_              = zeros(3, n);
            obj.omega_par_          = zeros(3, n);
            obj.int_dr_X_pv_r_      = zeros(n, 3);
            obj.int_pv_r_O_dd_r_    = zeros(n, 1);
            obj.grad_int_dr_        = zeros(n, 3);
            obj.grad_int_r_X_dr_    = zeros(n, 3);
            obj.grad_J_             = zeros(3, 3, n);
            obj.grad_v_com_         = zeros(n, 3);
            obj.K_                  = zeros(n, 1);
            obj.D_                  = zeros(n, 1);
        end
    end

    methods (Access = public)
        %Convert body to a struct
        function s = toStruct(obj)
            s = struct('BodyType', class(obj), 'BodyParameters', {{obj.Parameters}}, 'BodyDoF', obj.n);
        end
        
        %Update the body to the current configuration. 
        function Update(obj, q, dq, ddq)
            %Kinematic terms
            obj.T_                  = obj.T(q);
            obj.v_rel_              = obj.v_rel(q, dq);
            obj.omega_rel_          = obj.omega_rel(q, dq);
            obj.a_rel_              = obj.a_rel(q, dq, ddq);
            obj.domega_rel_         = obj.domega_rel(q, dq, ddq);
            obj.v_par_              = obj.v_par(q);
            obj.omega_par_          = obj.omega_par(q);
            %Inertial quantities
            obj.p_com_              = obj.p_com(q);
            obj.v_com_rel_          = obj.v_com_rel(q, dq);
            obj.a_com_rel_          = obj.a_com_rel(q, dq, ddq);
            obj.I_                  = obj.I(q);
            obj.m_                  = obj.m();
            obj.J_                  = obj.J(q, dq);
            obj.int_dr_             = obj.int_dr(q, dq);
            obj.int_ddr_            = obj.int_ddr(q, dq, ddq);
            obj.int_r_X_dr_         = obj.int_r_X_dr(q, dq);
            obj.int_r_X_ddr_        = obj.int_r_X_ddr(q, dq, ddq);
            obj.int_dr_X_pv_r_      = obj.int_dr_X_pv_r(q, dq);
            obj.int_pv_r_O_dd_r_    = obj.int_pv_r_O_dd_r(q, dq, ddq);
            obj.int_dr_O_dr_        = obj.int_dr_O_dr(q, dq);
            obj.grad_int_dr_        = obj.grad_int_dr(q);
            obj.grad_int_r_X_dr_    = obj.grad_int_r_X_dr(q);
            obj.grad_J_             = obj.grad_J(q);
            obj.grad_v_com_         = obj.grad_v_com(q);
            obj.K_                  = obj.K(q);
            obj.D_                  = obj.D(q, dq);
        end
    end
end



