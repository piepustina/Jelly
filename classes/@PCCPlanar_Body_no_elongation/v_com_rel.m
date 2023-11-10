function v_i_comi = v_com_rel(theta,dtheta,in3)
%V_COM_REL
%    V_I_COMI = V_COM_REL(THETA,DTHETA,IN3)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    06-Nov-2023 16:39:46

L_0 = in3(1,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = 1.0./theta.^2;
t5 = 1.0./theta.^3;
t6 = t2-1.0;
v_i_comi = [dtheta.*(L_0.*t4.*t6-L_0.*t5.*(t3-theta).*2.0);0.0;-dtheta.*(L_0.*t3.*t4+L_0.*t5.*t6.*2.0)];
end
