function out = J_r_i_X_dr_i_s(q,params,s)
out = skew(PACCylindrical_Delta_Body.r_i_s(q,params,s))*PACCylindrical_Delta_Body.J_dr_i_s(q,params,s)*PACCylindrical_Delta_Body.rho_L_s(params,s);
end
