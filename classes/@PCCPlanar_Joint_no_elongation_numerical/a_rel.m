function dv_i_1_i = a_rel(theta,dtheta,ddtheta,L_0)
%A_REL
%    DV_I_1_I = A_REL(THETA,DTHETA,DDTHETA,L_0)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    13-Oct-2023 18:22:41

t2 = theta./8.0;
t3 = theta./1.6e+1;
t4 = cos(t2);
t5 = cos(t3);
t6 = sin(t2);
t7 = sin(t3);
t8 = t5.^2;
t11 = t6.^2;
t9 = t8.^2;
t10 = t8.^3;
t12 = t11.^2;
t13 = t11.^3;
t15 = t8.^5;
t17 = t8.^7;
t19 = t11.*1.26e+2;
t21 = t8.*8.3e+1;
t14 = t9.^2;
t16 = t9.^3;
t18 = t13.*1.12e+2;
t20 = t12.*2.16e+2;
t23 = t9.*9.88e+2;
t24 = t10.*5.132e+3;
t26 = t17.*3.84e+3;
t29 = t15.*1.9072e+4;
t22 = -t20;
t27 = t14.*1.3568e+4;
t28 = t16.*1.3568e+4;
t30 = -t27;
t31 = -t28;
t32 = t18+t19+t22-2.1e+1;
dv_i_1_i = [dtheta.*((L_0.*dtheta.*t5.*(t5.*t7.*(8.3e+1./8.0)-t5.^3.*t7.*2.47e+2+t5.^5.*t7.*1.9245e+3-t5.^7.*t7.*6.784e+3+t5.^9.*t7.*1.192e+4-t5.^11.*t7.*1.0176e+4+t5.^13.*t7.*3.36e+3))./2.0+(L_0.*dtheta.*t7.*(t21-t23+t24+t26+t29+t30+t31-2.0))./3.2e+1)-(L_0.*ddtheta.*t5.*(t21-t23+t24+t26+t29+t30+t31-2.0))./2.0;0.0;dtheta.*((L_0.*dtheta.*t6.*(t4.*t6.*(6.3e+1./2.0)-t4.*t6.^3.*1.08e+2+t4.*t6.^5.*8.4e+1))./8.0+(L_0.*dtheta.*t4.*t32)./6.4e+1)+(L_0.*ddtheta.*t6.*t32)./8.0];
end
