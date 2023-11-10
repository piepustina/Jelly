function out = grad_J_s(q,params,s)
out = (sympagemtimes(pageskew(PCSCylindrical_Body.J_dr_i_s(q,params,s),1),skew(PCSCylindrical_Body.r_i_s(q,params,s)))+sympagemtimes(skew(PCSCylindrical_Body.r_i_s(q,params,s))',pageskew(PCSCylindrical_Body.J_dr_i_s(q,params,s),0)))*PCSCylindrical_Body.rho_L_s(params,s);
end
