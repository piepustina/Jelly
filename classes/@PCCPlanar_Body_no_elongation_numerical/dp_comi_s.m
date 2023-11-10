function out = dp_comi_s(q,dq,params,s)
out = (1/PCCPlanar_Body_no_elongation_numerical.m(params))*(PCCPlanar_Body_no_elongation_numerical.dp_i_s(q,dq,params,s))*(PCCPlanar_Body_no_elongation_numerical.rho_L_s(params,s));
end
