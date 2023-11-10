function out = v_com_rel(q,dq,params)
out = myintegral(@(s)PCSCylindrical_Body.dp_comi_s(q,dq,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
