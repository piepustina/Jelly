function out = ddp_comi_s(q,dq,ddq,params,s)
out = (1/PACCylindrical_Delta_Body.m(params))*(PACCylindrical_Delta_Body.ddp_i_s(q,dq,ddq,params,s))*(PACCylindrical_Delta_Body.rho_L_s(params,s));
end
