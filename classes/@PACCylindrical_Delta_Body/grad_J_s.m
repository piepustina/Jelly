function out = grad_J_s(q,params,s)
out = (sympagemtimes(pageskew(PACCylindrical_Delta_Body.J_dr_i_s(q,params,s),1),skew(PACCylindrical_Delta_Body.r_i_s(q,params,s)))+sympagemtimes(skew(PACCylindrical_Delta_Body.r_i_s(q,params,s))',pageskew(PACCylindrical_Delta_Body.J_dr_i_s(q,params,s),0)))*PACCylindrical_Delta_Body.rho_L_s(params,s);
end
