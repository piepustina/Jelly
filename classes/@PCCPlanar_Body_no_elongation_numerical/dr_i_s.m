function out = dr_i_s(q,dq,params,s)
out = PCCPlanar_Body_no_elongation_numerical.dp_i_s(q,dq,params,s)-myintegral(@(s)PCCPlanar_Body_no_elongation_numerical.dp_comi_s(q,dq,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
