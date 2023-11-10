function out = grad_J_s(q,params,s)
out = (sympagemtimes(pageskew(PCCPlanar_Body_no_elongation_numerical.J_dr_i_s(q,params,s),1),skew(PCCPlanar_Body_no_elongation_numerical.r_i_s(q,params,s)))+sympagemtimes(skew(PCCPlanar_Body_no_elongation_numerical.r_i_s(q,params,s))',pageskew(PCCPlanar_Body_no_elongation_numerical.J_dr_i_s(q,params,s),0)))*PCCPlanar_Body_no_elongation_numerical.rho_L_s(params,s);
end
