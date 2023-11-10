classdef PCSCylindrical_Joint < Joint
    %PCSCYLINDRICAL_JOINT Defines a joint for a isotropic PCS body
    
    properties(Constant)
        Type  = "PCSCylindrical";
        n = 6;
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
        function obj = PCSCylindrical_Joint(L_0)
            obj = obj@Joint(PCSCylindrical_Joint.n);
            obj.Parameters = L_0;
        end
    end
end

