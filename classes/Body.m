classdef Body < handle
    %Abstract class modeling a generic body of a :class:`BodyTree`. 

    properties (Access = public)
        T_                  = eye(4)
        v_rel_              = zeros(3, 1)
        omega_rel_          = zeros(3, 1)
        a_rel_              = zeros(3, 1)
        domega_rel_         = zeros(3, 1)
        v_par_
        omega_par_
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
        %Number of DoF of the body.
        n;
    end

    properties(Abstract)
        %Vector of parameters for the body.
        Parameters;
    end

    % Methods required to update the body state
    methods(Access = public)
            %% Kinematic terms
            function T_ = T(obj, q)
                %Evaluate the relative transformation matrix from the body distal end frame to the body proximal end frame.
                %
                %Args:
                %    q ([double], [sym]): Configuration variables

                T_ = eye(4, 'like', q);
            end

            function v_rel_ = v_rel(obj, q, dq)
                %Evaluate the relative linear velocity of the body distal end frame in the body proximal end frame.
                %
                %Args:
                %    q  ([double], [sym]): Configuration variables
                %    dq ([double], [sym]): First-order time derivative of the configuration variables

                v_rel_ = zeros(3, 1, 'like', q);
            end

            function omega_rel_ = omega_rel(obj, q, dq)
                %Evaluate the relative angular velocity of the body distal end frame in the body proximal end frame.
                %
                %Args:
                %    q  ([double], [sym]): Configuration variables
                %    dq ([double], [sym]): First-order time derivative of the configuration variables
                omega_rel_ = zeros(3, 1, 'like', q);
            end

            function a_rel_ = a_rel(obj, q, dq, ddq)
                %Evaluate the relative linear acceleration of the body distal end frame in the body proximal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of the configuration variables
                a_rel_ = zeros(3, 1, 'like', q);
            end
 
            function domega_rel_ = domega_rel(obj, q, dq, ddq)
                %Evaluate the relative angular acceleration of the body distal end frame in the body proximal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of the configuration variables
                domega_rel_ = zeros(3, 1, 'like', q);
            end

            function v_par_ = v_par(obj, q)
                %Evaluate the Jacobian of the relative linear velocity of the body distal end frame in its frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                v_par_ = zeros(3, obj.n, 'like', q);
            end
            
            function omega_par_ = omega_par(obj, q)
                %Evaluate the Jacobian of the relative angular velocity of the body distal end frame in its frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                omega_par_ = zeros(3, obj.n, 'like', q);
            end
            
            %% Inertial terms
            function p_com_ = p_com(obj, q)
                %Evaluate the center of mass position in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
            
                p_com_ = zeros(3, 1, 'like', q);
            end

            function v_com_rel_ = v_com_rel(obj, q, dq)
                %Evaluate the first-order time derivative of the center of mass position in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                v_com_rel_ = zeros(3, 1, 'like', q);
            end

            function a_com_rel_ = a_com_rel(obj, q, dq, ddq)
                %Evaluate the second-order time derivative of the center of mass position in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of the configuration variables
                a_com_rel_ = zeros(3, 1, 'like', q);
            end

            function I_ = I(obj, q)
                %Evaluate the body inertia matrix in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                I_ = zeros(3, 3, 'like', q);
            end

            function m_ = m(obj)
                %Evaluate the body mass.
                m_ = 0;
            end

            function J_ = J(obj, q, dq)
                %Evaluate the firs-order time derivative of the body inertia matrix in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                J_ = zeros(3, 3, 'like', q);
            end
            
            % Integral of \dot{r}
            % TODO: Always zero, remove.
            function int_dr_ = int_dr(~, q, ~)
                int_dr_ = zeros(3, 1, 'like', q);
            end
            % Integral of \ddot{r}
            % TODO: Always zero, remove.
            function int_ddr_ = int_ddr(~, q, ~, ~)
                int_ddr_ = zeros(3, 1, 'like', q);
            end

            % Integral of \cross(r, \dot{r})
            function int_r_X_dr_ = int_r_X_dr(obj, q, dq)
                %Evaluate :math:`\int_{V} r \times \dot{r} \,\, \mathrm{d}V`, where :math:`r` is the relative position vector of the body particle in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables

                int_r_X_dr_ = zeros(3, 1, 'like', q);
            end

            % Integral of \cross(r, \ddot{r})
            function int_r_X_ddr_ = int_r_X_ddr(obj, q, dq, ddq)
                %Evaluate :math:`\int_{V} r \times \ddot{r} \,\, \mathrm{d}V`, where :math:`r` is the relative position vector of the body particle in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of the configuration variables
                int_r_X_ddr_ = zeros(3, 1, 'like', q);
            end

            % Integral of \cross(\dor{r}, \jacobian{r}{q})
            function int_dr_X_pv_r_ = int_dr_X_pv_r(obj, q, dq)
                %Evaluate :math:`\int_{V} \left(\tilde{\dot{r}} \frac{\partial r}{\partial q}\right)^{T} \,\, \mathrm{d}V`, where :math:`r` is the relative position vector of the body particle in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                int_dr_X_pv_r_ = zeros(obj.n, 3, 'like', q);
            end

            % Integral of \dot(\jacobian{r}{q}, \ddot{r})
            function int_pv_r_O_dd_r_ = int_pv_r_O_dd_r(obj, q, dq, ddq)
                %Evaluate :math:`\int_{V} \left(\frac{\partial r}{\partial q}\right)^{T} \ddot{r} \,\, \mathrm{d}V`, where :math:`r` is the relative position vector of the body particle in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of the configuration variables
                int_pv_r_O_dd_r_ = zeros(obj.n, 1, 'like', q);
            end
            
            % Integral of \dot(\dot{r}, \dot{r})
            % TODO: Always zero, remove.
            function int_dr_O_dr_ = int_dr_O_dr(~, ~, ~)
                int_dr_O_dr_ = 0;
            end
            
            % Jacobian of the integral of \dot{r}
            % TODO: Always zero, remove.
            function grad_int_dr_ = grad_int_dr(obj, q)
                grad_int_dr_ = zeros(obj.n, 3, 'like', q);
            end

            % Jacobian of the integral of \cross{r, \dot{r}}
            function grad_int_r_X_dr_ = grad_int_r_X_dr(obj, q)
                %Evaluate :math:`\frac{\partial}{\partial q}\left(\int_{V} r \times \ddot{r} \,\, \mathrm{d}V\right)^{T}`, where :math:`r` is the relative position vector of the body particle in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                grad_int_r_X_dr_ = zeros(obj.n, 3, 'like', q);
            end

            % Jacobian of the time derivative of the inertia
            function grad_J_ = grad_J(obj, q)
                %Evaluate :math:`\frac{\partial \dot{J}}{\partial \dot{q}}`, where :math:`J` is the first-order time derivative of the body inertia.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                grad_J_ = zeros(3, 3, obj.n, 'like', q);
            end

            % Jacobian of the center of mass velocity
            function grad_v_com_ = grad_v_com(obj, q)
                %Evaluate :math:`\frac{\partial p_{\mathrm{CoM}}}{\partial q}`, where :math:`p_{\mathrm{CoM}}` is the position of the center of mass in the body distal end frame.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %
                grad_v_com_ = zeros(obj.n, 3, 'like', q);
            end

            % Generalized elastic force
            function K_ = K(obj, q)
                %Evaluate the generalized elastic force.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %
                K_ = zeros(obj.n, 1, 'like', q);
            end
            % Generalized damping force
            function D_ = D(obj, q, dq)
                %Evaluate the generalized damping force.
                %
                %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                D_ = zeros(obj.n, 1, 'like', q);
            end
    end

    methods (Access = protected)
        %Class constructor
        function obj = Body(n)
            %Construct the body.
            %
            %Args:
            %   n (double): Number of DoF of the body
            
            % Initialize all the variables
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
        function s = toStruct(obj)
            %Convert the body object to a struct representation.
            s = struct('BodyType', class(obj), 'BodyParameters', obj.Parameters, 'BodyDoF', obj.n);
        end
         
        function Update(obj, q, dq, ddq)
            %Update the current state of the body.
            %
            %Args:
                %    q   ([double], [sym]): Configuration variables
                %    dq  ([double], [sym]): First-order time derivative of the configuration variables
                %    ddq ([double], [sym]): Second-order time derivative of the configuration variables

            % Kinematic terms
            obj.T_                  = obj.T(q);
            obj.v_rel_              = obj.v_rel(q, dq);
            obj.omega_rel_          = obj.omega_rel(q, dq);
            obj.a_rel_              = obj.a_rel(q, dq, ddq);
            obj.domega_rel_         = obj.domega_rel(q, dq, ddq);
            obj.v_par_              = obj.v_par(q);
            obj.omega_par_          = obj.omega_par(q);
            % Inertial quantities
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



