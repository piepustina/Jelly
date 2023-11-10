function out = dr_i_s(q,dq,params,s)
out = PACCylindrical_Delta_Body.dp_i_s(q,dq,params,s)-myintegral(@(s)PACCylindrical_Delta_Body.dp_comi_s(q,dq,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
