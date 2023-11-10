function out = int_r_X_dr(q,dq,params)
out = myintegral(@(s)cross(PCCPlanar_Body_no_elongation_numerical.r_i_s(q,params,s),PCCPlanar_Body_no_elongation_numerical.dr_i_s(q,dq,params,s))*PCCPlanar_Body_no_elongation_numerical.rho_L_s(params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
