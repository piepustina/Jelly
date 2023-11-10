function out = J_r_i_X_dr_i_s(q,params,s)
out = skew(PCCPlanar_Body_no_elongation_numerical.r_i_s(q,params,s))*PCCPlanar_Body_no_elongation_numerical.J_dr_i_s(q,params,s)*PCCPlanar_Body_no_elongation_numerical.rho_L_s(params,s);
end
