function T = T(theta,in2)
%T
%    T = T(THETA,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-Jul-2023 11:50:46

a_ = in2(2,:);
alpha_ = in2(1,:);
d_ = in2(3,:);
t2 = cos(alpha_);
t3 = sin(alpha_);
t4 = cos(theta);
t5 = sin(theta);
T = reshape([t4,t5,0.0,0.0,-t2.*t5,t2.*t4,t3,0.0,t3.*t5,-t3.*t4,t2,0.0,a_.*t4,a_.*t5,d_,1.0],[4,4]);
end
