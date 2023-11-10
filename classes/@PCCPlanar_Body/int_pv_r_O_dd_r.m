function int_pv_r_O_dd_r = int_pv_r_O_dd_r(in1,in2,in3,in4)
%int_pv_r_O_dd_r
%    int_pv_r_O_dd_r = int_pv_r_O_dd_r(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-Jul-2023 11:47:47

L_0 = in4(1,:);
R = in4(2,:);
dddeltaL = in3(2,:);
ddeltaL = in2(2,:);
ddtheta = in3(1,:);
deltaL = in1(2,:);
dtheta = in2(1,:);
rho = in4(3,:);
theta = in1(1,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = R.^2;
t5 = dtheta.^2;
t6 = theta.^2;
t7 = theta.^3;
t9 = theta.^5;
t8 = t6.^2;
t10 = dddeltaL.*t8.*3.0;
t12 = L_0.*t5.*t8;
t13 = deltaL.*t5.*t8;
t11 = -t10;
mt1 = [(L_0.*rho.*t4.*1.0./theta.^7.*pi.*(L_0+deltaL).*(t11-t12-t13+L_0.*t5.*7.2e+1+dddeltaL.*t6.*1.2e+1+deltaL.*t5.*7.2e+1+L_0.*ddtheta.*t9-L_0.*ddtheta.*theta.*2.4e+1-L_0.*t2.*t5.*7.2e+1+ddtheta.*deltaL.*t9+ddeltaL.*dtheta.*t9.*2.0-ddtheta.*deltaL.*theta.*2.4e+1-ddeltaL.*dtheta.*theta.*4.8e+1-dddeltaL.*t2.*t6.*1.2e+1-dddeltaL.*t3.*t7.*3.0-deltaL.*t2.*t5.*7.2e+1+L_0.*t2.*t5.*t6.*6.0-L_0.*t3.*t5.*theta.*4.2e+1+ddtheta.*deltaL.*t3.*t6.*1.2e+1+ddeltaL.*dtheta.*t3.*t6.*2.4e+1+ddtheta.*deltaL.*t2.*theta.*2.4e+1+ddeltaL.*dtheta.*t2.*theta.*4.8e+1+deltaL.*t2.*t5.*t6.*6.0-deltaL.*t3.*t5.*theta.*4.2e+1+L_0.*ddtheta.*t3.*t6.*1.2e+1+L_0.*ddtheta.*t2.*theta.*2.4e+1))./3.0];
mt2 = [L_0.*rho.*t4.*1.0./t6.^3.*pi.*(t11+t12+t13+L_0.*t5.*3.6e+1+dddeltaL.*t6.*6.0+deltaL.*t5.*3.6e+1+L_0.*ddtheta.*t7.*3.0-L_0.*ddtheta.*theta.*1.2e+1-L_0.*t2.*t5.*3.6e+1-L_0.*t5.*t6.*9.0+ddtheta.*deltaL.*t7.*3.0+ddeltaL.*dtheta.*t7.*6.0-ddtheta.*deltaL.*theta.*1.2e+1-ddeltaL.*dtheta.*theta.*2.4e+1-dddeltaL.*t2.*t6.*6.0-deltaL.*t2.*t5.*3.6e+1-deltaL.*t5.*t6.*9.0+L_0.*t2.*t5.*t6.*3.0-L_0.*t3.*t5.*theta.*1.2e+1+ddtheta.*deltaL.*t3.*t6.*3.0+ddeltaL.*dtheta.*t3.*t6.*6.0+ddtheta.*deltaL.*t2.*theta.*1.2e+1+ddeltaL.*dtheta.*t2.*theta.*2.4e+1+deltaL.*t2.*t5.*t6.*3.0-deltaL.*t3.*t5.*theta.*1.2e+1+L_0.*ddtheta.*t3.*t6.*3.0+L_0.*ddtheta.*t2.*theta.*1.2e+1).*(-1.0./3.0)];
int_pv_r_O_dd_r = [mt1;mt2];
end
