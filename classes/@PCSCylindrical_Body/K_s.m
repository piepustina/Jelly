function K_s = K_s(in1,in2,s)
%K_s
%    K_s = K_s(IN1,IN2,S)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-Jul-2023 11:59:27

E_mod = in2(4,:);
L_0 = in2(1,:);
Poi = in2(5,:);
R = in2(2,:);
deltaL = in1(6,:);
gamma_X = in1(4,:);
gamma_Y = in1(5,:);
kappa_X = in1(1,:);
kappa_Y = in1(2,:);
tau_Z = in1(3,:);
t2 = Poi.*2.0;
t3 = R.^2;
t5 = 1.0./L_0;
t4 = t3.^2;
t6 = t5.^2;
t7 = t2+2.0;
t8 = 1.0./t7;
K_s = [(E_mod.*kappa_X.*t4.*t6.*pi)./4.0;(E_mod.*kappa_Y.*t4.*t6.*pi)./4.0;(E_mod.*t4.*t6.*t8.*tau_Z.*pi)./2.0;E_mod.*gamma_X.*t3.*t6.*t8.*pi;E_mod.*gamma_Y.*t3.*t6.*t8.*pi;E_mod.*t3.*t5.*pi.*(t5.*(L_0+deltaL)-1.0)];
end
