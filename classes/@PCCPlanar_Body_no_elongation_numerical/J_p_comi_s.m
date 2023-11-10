function J_p_comi_s = J_p_comi_s(theta,in2,s)
%J_p_comi_s
%    J_p_comi_s = J_p_comi_s(THETA,IN2,S)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    13-Oct-2023 18:22:37

L_0 = in2(1,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = s.^2;
t5 = 1.0./L_0;
t6 = theta./2.0;
t7 = theta./4.0;
t12 = theta./8.0;
t13 = theta./1.6e+1;
t8 = cos(t6);
t9 = cos(t7);
t10 = sin(t6);
t11 = sin(t7);
t14 = cos(t12);
t15 = cos(t13);
t16 = sin(t12);
t17 = sin(t13);
t20 = s.*t5.*t6;
t21 = s.*t5.*t7;
t26 = s.*t5.*t12;
t27 = s.*t5.*t13;
t18 = t14.^2;
t19 = t2.*t15;
t22 = cos(t20);
t23 = cos(t21);
t24 = sin(t20);
t25 = sin(t21);
t29 = cos(t26);
t30 = cos(t27);
t31 = sin(t26);
t32 = sin(t27);
t28 = t18.*2.0;
t33 = -t19;
t35 = t29.^2;
t34 = t28-1.0;
t36 = t35.*2.0;
t37 = t2+t33+1.0;
t38 = t34.^2;
t40 = t36-1.0;
t39 = t38.*2.0;
t42 = t40.^2;
t41 = t39-1.0;
t43 = t42.*2.0;
t44 = t43-1.0;
mt1 = [t5.*(L_0.*t8.*t9.*t14.*t37.*(-1.0./2.0)+(L_0.*t9.*t10.*t16.*t37)./8.0+(L_0.*t10.*t11.*t14.*t37)./4.0+s.*t2.*t29.*t40.*t44-L_0.*t9.*t10.*t14.*(-t3+t3.*t15+(t2.*t17)./1.6e+1)+s.*t3.*t23.*t24.*t29.*t30-t3.*t4.*t5.*t31.*t35.*t42.*2.0-(t3.*t4.*t5.*t31.*t35.*t44)./2.0-(t3.*t4.*t5.*t31.*t40.*t44)./8.0-(t2.*t4.*t5.*t22.*t23.*t29.*t30)./2.0+(t2.*t4.*t5.*t23.*t24.*t29.*t32)./1.6e+1+(t2.*t4.*t5.*t23.*t24.*t30.*t31)./8.0+(t2.*t4.*t5.*t24.*t25.*t29.*t30)./4.0);0.0];
mt2 = [t5.*(L_0.*t9.*t10.*t14.*t33+(L_0.*t2.*t16.*t18.*t41)./2.0+L_0.*t2.*t16.*t28.*t38+L_0.*t3.*t14.*t34.*t41+(L_0.*t2.*t16.*t34.*t41)./8.0-s.*t3.*t29.*t40.*t44+s.*t2.*t23.*t24.*t29.*t30-t2.*t4.*t5.*t31.*t35.*t42.*2.0-(t2.*t4.*t5.*t31.*t35.*t44)./2.0-(t2.*t4.*t5.*t31.*t40.*t44)./8.0-(L_0.*t3.*t8.*t9.*t14.*t15)./2.0+(L_0.*t3.*t9.*t10.*t14.*t17)./1.6e+1+(L_0.*t3.*t9.*t10.*t15.*t16)./8.0+(L_0.*t3.*t10.*t11.*t14.*t15)./4.0+(t3.*t4.*t5.*t22.*t23.*t29.*t30)./2.0-(t3.*t4.*t5.*t23.*t24.*t29.*t32)./1.6e+1-(t3.*t4.*t5.*t23.*t24.*t30.*t31)./8.0-(t3.*t4.*t5.*t24.*t25.*t29.*t30)./4.0)];
J_p_comi_s = [mt1;mt2];
end
