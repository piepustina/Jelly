function I = I(in1,in2)
%I
%    I = I(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-Jul-2023 11:47:46

L_0 = in2(1,:);
R = in2(2,:);
deltaL = in1(2,:);
rho = in2(3,:);
theta = in1(1,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = L_0+deltaL;
t5 = R.^2;
t6 = theta.*2.0;
t7 = theta.^2;
t8 = t2.*2.0;
t9 = sin(t6);
t10 = t3.*theta;
t11 = 1.0./t7.^2;
t12 = t4.^2;
t13 = t8+t10-2.0;
t14 = (L_0.*rho.*t3.*t5.*t11.*t12.*t13.*pi)./2.0;
t15 = -t14;
I = reshape([(L_0.*rho.*t5.*t11.*t12.*pi.*(cos(t6).*2.0+t6.*theta+t9.*theta-2.0))./4.0,t15,0.0,t15,L_0.*rho.*t5.*t11.*t12.*pi.*(t2.*-4.0-t7+t2.*t8+(t9.*theta)./2.0+2.0).*(-1.0./2.0),0.0,0.0,0.0,L_0.*rho.*t5.*t11.*t12.*pi.*(t7+t8-2.0)],[3,3]);
end
