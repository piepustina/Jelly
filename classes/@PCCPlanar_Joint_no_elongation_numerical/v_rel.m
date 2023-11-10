function v_i_1_i = v_rel(theta,dtheta,L_0)
%V_REL
%    V_I_1_I = V_REL(THETA,DTHETA,L_0)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    13-Oct-2023 18:22:41

t2 = theta./8.0;
t3 = theta./1.6e+1;
t4 = cos(t3);
t5 = sin(t2);
v_i_1_i = [L_0.*dtheta.*t4.*(t4.^2.*8.3e+1-t4.^4.*9.88e+2+t4.^6.*5.132e+3-t4.^8.*1.3568e+4+t4.^10.*1.9072e+4-t4.^12.*1.3568e+4+t4.^14.*3.84e+3-2.0).*(-1.0./2.0);0.0;(L_0.*dtheta.*t5.*(t5.^2.*1.26e+2-t5.^4.*2.16e+2+t5.^6.*1.12e+2-2.1e+1))./8.0];
end
