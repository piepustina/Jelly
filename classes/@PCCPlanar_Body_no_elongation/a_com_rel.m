function dv_i_comi = a_com_rel(theta,dtheta,ddtheta,in4)
%A_COM_REL
%    DV_I_COMI = A_COM_REL(THETA,DTHETA,DDTHETA,IN4)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    06-Nov-2023 16:39:46

L_0 = in4(1,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = dtheta.^2;
t5 = 1.0./theta.^2;
t6 = 1.0./theta.^3;
t7 = t5.^2;
t8 = t2-1.0;
t11 = L_0.*t3.*t5;
dv_i_comi = [-t4.*(t11+L_0.*t6.*t8.*4.0-L_0.*t7.*(t3-theta).*6.0)+ddtheta.*(L_0.*t5.*t8-L_0.*t6.*(t3-theta).*2.0);0.0;-ddtheta.*(t11+L_0.*t6.*t8.*2.0)+t4.*(-L_0.*t2.*t5+L_0.*t3.*t6.*4.0+L_0.*t7.*t8.*6.0)];
end
