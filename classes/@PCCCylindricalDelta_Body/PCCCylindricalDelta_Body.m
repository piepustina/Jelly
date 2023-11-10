classdef PCCCylindricalDelta_Body < Body
    %PCCCYLINDRICAL_BODY Defines a linear isotropic PCC body with elongation.
    
    properties(Constant)
        Type = "PCCCylindricalDelta";
        n    = 3;
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
    end

    methods
        function obj = PCCCylindricalDelta_Body(Parameters)
            if isrow(Parameters)
                Parameters = Parameters';
            end
            obj = obj@Body(PCCCylindricalDelta_Body.n);
            obj.Parameters = Parameters;
            obj.RestLength = Parameters(1);
        end

        
    end
end

