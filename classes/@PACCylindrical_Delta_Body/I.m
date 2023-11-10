function out = I(q,params)
out = myintegral(@(s)skew(PACCylindrical_Delta_Body.r_i_s(q,params,s))'*skew(PACCylindrical_Delta_Body.r_i_s(q,params,s))*PACCylindrical_Delta_Body.rho_L_s(params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
