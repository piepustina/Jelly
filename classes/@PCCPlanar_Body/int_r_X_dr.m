function int_r_i_X_dr_i = int_r_X_dr(in1,in2,in3)
%int_r_X_dr
%    int_r_i_X_dr_i = int_r_X_dr(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-Jul-2023 11:47:46

L_0 = in3(1,:);
R = in3(2,:);
deltaL = in1(2,:);
dtheta = in2(1,:);
rho = in3(3,:);
theta = in1(1,:);
int_r_i_X_dr_i = [0.0;0.0;L_0.*R.^2.*dtheta.*rho.*1.0./theta.^4.*pi.*(L_0+deltaL).^2.*(cos(theta).*2.0+theta.^2-2.0).*(-1.0./2.0)];
end
