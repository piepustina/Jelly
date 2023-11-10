classdef JointOld < handle
    %JOINT Class that represents a generic Joint
    %#codegen

    %TODO: A joint should be an element of a Body since each body is
    %attached through a joint to the remaining part of the system.
    
    properties(Access = public)
        T_          = eye(4)
        v_rel_      = zeros(3, 1)
        omega_rel_  = zeros(3, 1)
        a_rel_      = zeros(3, 1)
        domega_rel_ = zeros(3, 1)
        v_par_
        omega_par_
    end

    properties (Abstract, Constant)
        %Joint type
        Type
        %Number of DOFs allowed by the joint
        n
    end

    properties(Abstract)
        %Vector containing the parameters of the joint
        Parameters
    end

    %Every class that inherits has to implement these methods for the
    %update.
    methods (Abstract, Static)
            T_          = T(q, Parameters);
            T_s_        = T_s(q, Parameters, s);
            v_rel_      = v_rel(q, dq, Parameters);
            omega_rel_  = omega_rel(q, dq, Parameters);
            a_rel_      = a_rel(q, dq, ddq, Parameters);
            domega_rel_ = domega_rel(q, dq, ddq, Parameters);
            v_par_      = v_par(q, Parameters);
            omega_par_  = omega_par(q, Parameters);
    end

    methods (Access = protected)
        function obj = JointOld(n)
            obj.v_par_     = zeros(3, n);
            obj.omega_par_ = zeros(3, n);
        end
    end

    methods (Access = public)
        %Provide a struct representation of the Joint
        function s = toStruct(obj)
            s = struct('JointType', class(obj), 'JointParameters', {{obj.Parameters}}, 'JointDoF', obj.n);
        end

        %Update the joint parameters into the current configuration
        function updateJoint(obj, q, dq, ddq)
            obj.T_          = obj.T(q, obj.Parameters);
            obj.v_rel_      = obj.v_rel(q, dq, obj.Parameters);
            obj.omega_rel_  = obj.omega_rel(q, dq, obj.Parameters);
            obj.a_rel_      = obj.a_rel(q, dq, ddq, obj.Parameters);
            obj.domega_rel_ = obj.domega_rel(q, dq, ddq, obj.Parameters);
            obj.v_par_      = obj.v_par(q, obj.Parameters);
            obj.omega_par_  = obj.omega_par(q, obj.Parameters);
        end
    end


end


