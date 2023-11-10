classdef BodyNew < handle
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
        %Type;
        %Number of DoFs
        n;
    end

    properties(Abstract)
        %Dynamic parameters
        Parameters;
        %Rest length
        RestLength;
    end

    %Methods required to update the body state
    methods(Abstract)
            %Kinematic terms
            T_                  = T(q);
            T_s_                = T_s(q, s);
            v_rel_              = v_rel(q, dq);
            omega_rel_          = omega_rel(q, dq);
            a_rel_              = a_rel(q, dq, ddq);
            domega_rel_         = domega_rel(q, dq, ddq);
            v_par_              = v_par(q);
            omega_par_          = omega_par(q);
            %Inertial terms
            p_com_              = p_com(q);
            v_com_rel_          = v_com_rel(q, dq);
            a_com_rel_          = a_com_rel(q, dq, ddq);
            I_                  = I(q);
            m_                  = m();
            J_                  = J(q, dq);
            int_dr_             = int_dr(q, dq);
            int_ddr_            = int_ddr(q, dq, ddq);
            int_r_X_dr_         = int_r_X_dr(q, dq);
            int_r_X_ddr_        = int_r_X_ddr(q, dq, ddq);
            int_dr_X_pv_r_      = int_dr_X_pv_r(q, dq);
            int_pv_r_O_dd_r_    = int_pv_r_O_dd_r(q, dq, ddq);
            int_dr_O_dr_        = int_dr_O_dr(q, dq);
            grad_int_dr_        = grad_int_dr(q);
            grad_int_r_X_dr_    = grad_int_r_X_dr(q);
            grad_J_             = grad_J(q);
            grad_v_com_         = grad_v_com(q);
            xi_                 = xi(q, s);
            %Generalized elastic and damping forces
            K_ = K(q);
            D_ = D(q, dq);
    end

    methods (Access = protected)
        function obj = BodyNew(n)
            %Initialize the terms
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
        function updateBody(obj, q, dq, ddq)
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



