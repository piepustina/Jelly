function grad_v_com_i = grad_v_com(in1,in2)
%GRAD_V_COM
%    GRAD_V_COM_I = GRAD_V_COM(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-Jul-2023 11:47:47

L_0 = in2(1,:);
deltaL = in1(2,:);
theta = in1(1,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = L_0+deltaL;
t5 = 1.0./theta.^2;
t6 = 1.0./theta.^3;
t7 = t2-1.0;
t8 = -t3;
grad_v_com_i = reshape([t4.*t5.*t8-t4.*t6.*t7.*2.0,t5.*t7,-t4.*t5.*t7+t4.*t6.*(t3-theta).*2.0,-t5.*(t3-theta),0.0,0.0],[2,3]);
end
