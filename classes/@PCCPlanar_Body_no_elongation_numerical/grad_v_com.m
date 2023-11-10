function out = grad_v_com(q,params)
out = myintegral(@(s)PCCPlanar_Body_no_elongation_numerical.J_p_comi_s(q,params,s)',0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
