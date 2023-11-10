function out = grad_int_r_X_dr(q,params)
out = myintegral(@(s)PCCPlanar_Body_no_elongation_numerical.J_r_i_X_dr_i_s(q,params,s)',0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
