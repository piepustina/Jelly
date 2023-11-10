function out = grad_J(q,params)
out = myintegral(@(s)PACCylindrical_Delta_Body.grad_J_s(q,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
