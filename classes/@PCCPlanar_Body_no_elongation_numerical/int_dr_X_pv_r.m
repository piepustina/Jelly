function out = int_dr_X_pv_r(q,dq,params)
out = myintegral(@(s)PCCPlanar_Body_no_elongation_numerical.J_dr_i_s(q,params,s)'*skew(PCCPlanar_Body_no_elongation_numerical.dr_i_s(q,dq,params,s))'*PCCPlanar_Body_no_elongation_numerical.rho_L_s(params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
