function out = ddp_comi_s(q,dq,ddq,params,s)
out = (1/PCCPlanar_Body_no_elongation_numerical.m(params))*(PCCPlanar_Body_no_elongation_numerical.ddp_i_s(q,dq,ddq,params,s))*(PCCPlanar_Body_no_elongation_numerical.rho_L_s(params,s));
end
