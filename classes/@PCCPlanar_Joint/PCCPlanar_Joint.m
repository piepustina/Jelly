classdef PCCPlanar_Joint < Joint
    %PCC_PLANAR_SLENDER_NO_ELONGATION_JOINT

    properties(Constant)
        Type  = "PCC_planar";
        n = 2;
    end

    properties
        Parameters;
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

    methods
        function obj = PCCPlanar_Joint(L_0)
            obj = obj@Joint(PCCPlanar_Joint.n);
            obj.Parameters = L_0;
        end
    end
end

