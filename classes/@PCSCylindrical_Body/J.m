function out = J(q,dq,params)
out = myintegral(@(s)skew(PCSCylindrical_Body.dr_i_s(q,dq,params,s))'*skew(PCSCylindrical_Body.r_i_s(q,params,s))*PCSCylindrical_Body.rho_L_s(params,s)+skew(PCSCylindrical_Body.r_i_s(q,params,s))'*skew(PCSCylindrical_Body.dr_i_s(q,dq,params,s))*PCSCylindrical_Body.rho_L_s(params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
