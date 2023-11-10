function out = J_r_i_X_dr_i_s(q,params,s)
out = skew(PCSCylindrical_Body.r_i_s(q,params,s))*PCSCylindrical_Body.J_dr_i_s(q,params,s)*PCSCylindrical_Body.rho_L_s(params,s);
end
