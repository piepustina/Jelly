function xi_i = xi(theta,in2,s)
%XI
%    XI_I = XI(THETA,IN2,S)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    13-Oct-2023 18:12:57

L_0 = in2(1,:);
xi_i = [0.0;-theta./L_0;0.0;0.0;0.0;1.0];
end
