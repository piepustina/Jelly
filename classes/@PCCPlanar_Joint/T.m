function T = T(in1,L_0)
%T
%    T = T(IN1,L_0)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-Jul-2023 11:47:45

deltaL = in1(2,:);
theta = in1(1,:);
t2 = cos(theta);
t3 = sin(theta);
t4 = L_0+deltaL;
t5 = 1.0./theta;
T = reshape([t2,t3,0.0,0.0,-t3,t2,0.0,0.0,0.0,0.0,1.0,0.0,t3.*t4.*t5,-t4.*t5.*(t2-1.0),0.0,1.0],[4,4]);
end
