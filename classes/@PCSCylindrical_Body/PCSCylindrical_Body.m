classdef PCSCylindrical_Body < Body
    %PCSCYLINDRICAL_BODY Defines a linear isotropic PCS body.
    
    properties(Constant)
        Type = "PCSCylindrical";
        n    = 6;
    end

    properties
        %Parameters are in this order:
        %L_0 = rest length
        %R   = radius
        %rho = mass density
        %E   = Young's modulus
        %Poi = Poisson's ratio
        %Eta = Material damping
        Parameters;
        RestLength = 0;
    end

    %Implementation of the methods required by Body class
    methods(Static)
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
            K_               = K(q, Parameters);
            D_               = D(q, dq, Parameters);
            %% Functions to add when the Body is computed numerically
            r_i_s_ = r_i_s(q, Parameters, s);
            dr_i_s_ = dr_i_s(q, dq, Parameters, s);
            ddr_i_s_ = ddr_i_s(q, dq, ddq, Parameters, s);
            p_i_ = p_i_s(q, Parameters, s);
            dp_i_ = dp_i_s(q, dq, Parameters, s);
            ddp_i_ = ddp_i_s(q, dq, ddq, Parameters, s);
            p_comi_s_ = p_comi_s(q, Parameters, s);
            dp_comi_s_ = dp_comi_s(q, dq, Parameters, s);
            ddp_comi_s_ = ddp_comi_s(q, dq, ddq, Parameters, s);
            J_p_i_ = J_p_i(q, Parameters, s);
            J_p_comi_s_ = J_p_comi_s(q, Parameters, s);
            rho_L_s_ = rho_L_s(Parameters, s);
            J_dr_i_s_ = J_dr_i_s(q, Parameters, s);
            J_r_i_X_dr_i_s_ = J_r_i_X_dr_i_s(q, Parameters, s);
            grad_J_s_ = grad_J_s(q, Parameters, s);
            K_s_ = K_s(q, Parameters, s);
            D_s_ = D_s(q, dq, Parameters, s);
    end

    methods
        function obj = PCSCylindrical_Body(Parameters)
            if isrow(Parameters)
                Parameters = Parameters';
            end
            obj = obj@Body(PCSCylindrical_Body.n);
            obj.Parameters = Parameters;
            obj.RestLength = Parameters(1);
        end
    end
end

