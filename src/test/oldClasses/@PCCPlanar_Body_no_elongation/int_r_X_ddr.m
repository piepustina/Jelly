function int_r_i_X_ddr_i = int_r_X_ddr(theta,dtheta,ddtheta,in4)
%int_r_X_ddr
%    int_r_i_X_ddr_i = int_r_X_ddr(THETA,DTHETA,DDTHETA,IN4)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    06-Nov-2023 16:39:46

L_0 = in4(1,:);
R = in4(2,:);
rho = in4(3,:);
t2 = dtheta.^2;
t3 = theta./2.0;
t4 = sin(t3);
t5 = t4.^2;
int_r_i_X_ddr_i = [0.0;L_0.^3.*R.^2.*rho.*1.0./theta.^5.*pi.*(t2.*t5.*-1.6e+1-ddtheta.*theta.^3+t2.*theta.^2.*2.0+t2.*theta.*sin(theta).*2.0+ddtheta.*t5.*theta.*4.0).*(-1.0./2.0);0.0];
end