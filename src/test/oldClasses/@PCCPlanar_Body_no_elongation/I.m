function I = I(theta,in2)
%I
%    I = I(THETA,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    06-Nov-2023 16:39:46

L_0 = in2(1,:);
R = in2(2,:);
rho = in2(3,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = L_0.^3;
t5 = R.^2;
t6 = theta.*2.0;
t7 = theta.^2;
t8 = t2.*2.0;
t9 = sin(t6);
t10 = t3.*theta;
t11 = 1.0./t7.^2;
t12 = t8+t10-2.0;
t13 = (rho.*t3.*t4.*t5.*t11.*t12.*pi)./2.0;
I = reshape([rho.*t4.*t5.*t11.*pi.*(t2.*-4.0-t7+t2.*t8+(t9.*theta)./2.0+2.0).*(-1.0./2.0),0.0,t13,0.0,rho.*t4.*t5.*t11.*pi.*(t7+t8-2.0),0.0,t13,0.0,(rho.*t4.*t5.*t11.*pi.*(cos(t6).*2.0+t6.*theta+t9.*theta-2.0))./4.0],[3,3]);
end