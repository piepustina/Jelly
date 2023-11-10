function out = J_dr_i_s(q,params,s)
out = PACCylindrical_Delta_Body.J_p_i(q,params,s)-myintegral(@(s)PACCylindrical_Delta_Body.J_p_comi_s(q,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
