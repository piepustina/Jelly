function out = r_i_s(q,params,s)
out = PCSCylindrical_Body.p_i_s(q,params,s)-myintegral(@(s)PCSCylindrical_Body.p_comi_s(q,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
