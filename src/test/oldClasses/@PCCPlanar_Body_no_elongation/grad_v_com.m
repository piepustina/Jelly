function grad_v_com_i = grad_v_com(theta,in2)
%GRAD_V_COM
%    GRAD_V_COM_I = GRAD_V_COM(THETA,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    06-Nov-2023 16:39:47

L_0 = in2(1,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = 1.0./theta.^2;
t5 = 1.0./theta.^3;
t6 = t2-1.0;
grad_v_com_i = [L_0.*t4.*t6-L_0.*t5.*(t3-theta).*2.0,0.0,-L_0.*t3.*t4-L_0.*t5.*t6.*2.0];
end