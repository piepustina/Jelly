function v_i_comi = v_com_rel(in1,in2,in3)
%V_COM_REL
%    V_I_COMI = V_COM_REL(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-Jul-2023 11:47:46

L_0 = in3(1,:);
ddeltaL = in2(2,:);
deltaL = in1(2,:);
dtheta = in2(1,:);
theta = in1(1,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = L_0+deltaL;
t5 = 1.0./theta.^2;
t6 = 1.0./theta.^3;
t7 = t2-1.0;
v_i_comi = [-dtheta.*(t3.*t4.*t5+t4.*t6.*t7.*2.0)+ddeltaL.*t5.*t7;-dtheta.*(t4.*t5.*t7-t4.*t6.*(t3-theta).*2.0)-ddeltaL.*t5.*(t3-theta);0.0];
end
