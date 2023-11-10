function out = D(q,dq,params)
out = myintegral(@(s)PACCylindrical_Delta_Body.D_s(q,dq,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
