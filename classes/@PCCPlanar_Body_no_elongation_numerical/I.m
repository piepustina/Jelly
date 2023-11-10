function out = I(q,params)
out = myintegral(@(s)skew(PCCPlanar_Body_no_elongation_numerical.r_i_s(q,params,s))'*skew(PCCPlanar_Body_no_elongation_numerical.r_i_s(q,params,s))*PCCPlanar_Body_no_elongation_numerical.rho_L_s(params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
