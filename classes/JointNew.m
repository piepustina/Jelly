classdef JointNew < BodyNew
    %JOINT Class that represents a generic Joint. The class is represented
    %as a body with no inertial parameters.
    %#codegen
    
    % properties(Access = public)
    %     T_          = eye(4)
    %     v_rel_      = zeros(3, 1)
    %     omega_rel_  = zeros(3, 1)
    %     a_rel_      = zeros(3, 1)
    %     domega_rel_ = zeros(3, 1)
    %     v_par_
    %     omega_par_
    % end

    % properties (Abstract, Constant)
    %     %Number of degrees of freedom of the joint
    %     n
    % end

    % properties(Abstract)
    %     %Parameters of the joint
    %     Parameters
    % end

    % %Every class that inherits has to implement these methods for the
    % %update.
    % methods (Abstract, Static)
    %         T_          = T(q, Parameters);
    %         T_s_        = T_s(q, Parameters, s);
    %         v_rel_      = v_rel(q, dq, Parameters);
    %         omega_rel_  = omega_rel(q, dq, Parameters);
    %         a_rel_      = a_rel(q, dq, ddq, Parameters);
    %         domega_rel_ = domega_rel(q, dq, ddq, Parameters);
    %         v_par_      = v_par(q, Parameters);
    %         omega_par_  = omega_par(q, Parameters);
    % end

    methods (Access = protected)
        function obj = JointNew(n)
            obj = obj@BodyNew(n);
        end
    end

    methods (Access = public)
        %Overload the structure representation
        function s = toStruct(obj)
            s = struct('JointType', class(obj), 'JointParameters', {{obj.Parameters}}, 'JointDoF', obj.n);
        end
    end
end


