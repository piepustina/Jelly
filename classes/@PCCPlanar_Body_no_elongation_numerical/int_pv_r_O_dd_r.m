function out = int_pv_r_O_dd_r(q,dq,ddq,params)
out = myintegral(@(s)PCCPlanar_Body_no_elongation_numerical.J_dr_i_s(q,params,s)'*PCCPlanar_Body_no_elongation_numerical.ddr_i_s(q,dq,ddq,params,s)*PCCPlanar_Body_no_elongation_numerical.rho_L_s(params,s),0,params(1),'ArrayValued',1,'AbsTol',0.000100,'RelTol',0.001000);
end
