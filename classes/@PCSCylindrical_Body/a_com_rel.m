function out = a_com_rel(q,dq,ddq,params)
out = myintegral(@(s)PCSCylindrical_Body.ddp_comi_s(q,dq,ddq,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
