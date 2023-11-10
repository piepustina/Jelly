function out = a_com_rel(q,dq,ddq,params)
out = myintegral(@(s)PCCPlanar_Body_no_elongation_numerical.ddp_comi_s(q,dq,ddq,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
