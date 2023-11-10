function out = r_i_s(q,params,s)
out = PCCPlanar_Body_no_elongation_numerical.p_i_s(q,params,s)-myintegral(@(s)PCCPlanar_Body_no_elongation_numerical.p_comi_s(q,params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
