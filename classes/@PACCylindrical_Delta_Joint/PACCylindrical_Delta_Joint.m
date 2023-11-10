classdef PACCylindrical_Delta_Joint < Joint
    %PACCYLINDRICAL_DELTA_JOINT Defines a joint for a isotropic PAC body (thin rod) with
    %elongation
    
    properties(Constant)
        Type  = "PACCylindrical_Delta";
        n = 5;
    end

    %Methods implemented as separate files.
    methods (Static)
            T_          = T(q, Parameters);
            T_s         = T_s(q, Parameters, s);
            v_rel_      = v_rel(q, dq, Parameters);
            omega_rel_  = omega_rel(q, dq, Parameters);
            a_rel_      = a_rel(q, dq, ddq, Parameters);
            domega_rel_ = domega_rel(q, dq, ddq, Parameters);
            v_par_      = v_par(q, Parameters);
            omega_par_  = omega_par(q, Parameters);
    end

    properties
        %The parameters are:
        %- L_0 = rest length of the body
        Parameters;
    end

    methods
        function obj = PACCylindrical_Delta_Joint(L_0)
            obj = obj@Joint(PACCylindrical_Delta_Joint.n);
            obj.Parameters = L_0;
        end
    end
end

