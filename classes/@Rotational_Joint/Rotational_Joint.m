classdef Rotational_Joint < Joint
    %Rotational Joint

    properties(Constant)
        Type  = "Rotational_Joint___";
        n = 1;
    end

    properties
        Parameters;
        %Parameters are in this order:
        %alpha_
        %a_
        %d_
        %m_
        %pcom_x
        %pcom_y
        %pcom_z
        %I_x_x
        %I_y_y
        %I_z_z
        %I_x_y
        %I_x_z
        %I_y_z
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
        function obj = Rotational_Joint(Parameters)
            obj = obj@Joint(Rotational_Joint.n);
            obj.Parameters = Parameters;
        end
    end
end

