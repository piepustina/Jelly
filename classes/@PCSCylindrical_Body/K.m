function out = K(q,params)
out = myintegral(@(s)PCSCylindrical_Body.K_s(q,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
