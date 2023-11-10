classdef Body < handle
    %BODY 
    properties (Access = public)
        p_com_              = zeros(3, 1)
        v_com_rel_          = zeros(3, 1)
        a_com_rel_          = zeros(3, 1)
        m_                  = 0
        I_                  = zeros(3, 3)
        dI_                 = zeros(3, 3)
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
        Type;
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
    methods(Abstract, Static)
            p_com_           = p_com(q, Parameters);
            v_com_rel_       = v_com_rel(q, dq, Parameters);
            a_com_rel_       = a_com_rel(q, dq, ddq, Parameters);
            I_               = I(q, Parameters);
            m_               = m(Parameters);
            dI_              = dI(q, dq, Parameters);
            J_               = J(q, dq, Parameters);
            int_dr_          = int_dr(q, dq, Parameters);
            int_ddr_         = int_ddr(q, dq, ddq, Parameters);
            int_r_X_dr_      = int_r_X_dr(q, dq, Parameters);
            int_r_X_ddr_     = int_r_X_ddr(q, dq, ddq, Parameters);
            int_dr_X_pv_r_   = int_dr_X_pv_r(q, dq, Parameters);
            int_pv_r_O_dd_r_ = int_pv_r_O_dd_r(q, dq, ddq, Parameters);
            int_dr_O_dr_     = int_dr_O_dr(q, dq, Parameters);
            grad_int_dr_     = grad_int_dr(q, Parameters);
            grad_int_r_X_dr_ = grad_int_r_X_dr(q, Parameters);
            grad_J_          = grad_J(q, Parameters);
            grad_v_com_      = grad_v_com(q, Parameters);
            xi_              = xi(q, Parameters, s);
            %Generalized elastic and damping forces
            K_ = K(q, Parameters);
            D_ = D(q, dq, Parameters);
    end

    methods (Access = protected)
        function obj = Body(n)
            obj.int_dr_X_pv_r_      = zeros(n, 3);
            obj.int_pv_r_O_dd_r_    = zeros(n, 1);
            obj.grad_int_dr_        = zeros(n, 3);
            obj.grad_int_r_X_dr_    = zeros(n, 3);
            obj.grad_J_             = zeros(3, 3, n);
            obj.grad_v_com_         = zeros(n, 3);
            %Generalized elastic and damping forces
            obj.K_                  = zeros(n, 1);
            obj.D_                  = zeros(n, 1);
        end
    end

    methods (Access = public)
        %Convert body to a struct
        function s = toStruct(obj)
            s = struct('BodyType', class(obj), 'BodyParameters',  {{obj.Parameters}}, 'BodyDoF', obj.n);
        end
        
        %Update the body to the current configuration. 
        function updateBody(obj, q, dq, ddq)
            %Inertial quantities
            obj.p_com_           = obj.p_com(q, obj.Parameters);
            obj.v_com_rel_       = obj.v_com_rel(q, dq, obj.Parameters);
            obj.a_com_rel_       = obj.a_com_rel(q, dq, ddq, obj.Parameters);
            obj.I_               = obj.I(q, obj.Parameters);
            obj.m_               = obj.m(obj.Parameters);
            obj.dI_              = obj.dI(q, dq, obj.Parameters);
            obj.J_               = obj.J(q, dq, obj.Parameters);
            obj.int_dr_          = obj.int_dr(q, dq, obj.Parameters);
            obj.int_ddr_         = obj.int_ddr(q, dq, ddq, obj.Parameters);
            obj.int_r_X_dr_      = obj.int_r_X_dr(q, dq, obj.Parameters);
            obj.int_r_X_ddr_     = obj.int_r_X_ddr(q, dq, ddq, obj.Parameters);
            obj.int_dr_X_pv_r_   = obj.int_dr_X_pv_r(q, dq, obj.Parameters);
            obj.int_pv_r_O_dd_r_ = obj.int_pv_r_O_dd_r(q, dq, ddq, obj.Parameters);
            obj.int_dr_O_dr_     = obj.int_dr_O_dr(q, dq, obj.Parameters);
            obj.grad_int_dr_     = obj.grad_int_dr(q, obj.Parameters);
            obj.grad_int_r_X_dr_ = obj.grad_int_r_X_dr(q, obj.Parameters);
            obj.grad_J_          = obj.grad_J(q, obj.Parameters);
            obj.grad_v_com_      = obj.grad_v_com(q, obj.Parameters);
            obj.K_               = obj.K(q, obj.Parameters);
            obj.D_               = obj.D(q, dq, obj.Parameters);
        end
    end
end



