function out = ddr_i_s(q,dq,ddq,params,s)
out = PACCylindrical_Delta_Body.ddp_i_s(q,dq,ddq,params,s)-myintegral(@(s)PACCylindrical_Delta_Body.ddp_comi_s(q,dq,ddq,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
