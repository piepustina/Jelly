%% Export all the functions for the current model, those that are numerical separately
%--------------------------------------------------------------------------
% Store all the terms describing the body also in a common folder ---------
%--------------------------------------------------------------------------
disp("Exporting functions...");
cF = pwd;
cd(destination_folder_j);

matlabFunction(T                   , 'Vars', {q, [L_0]}, 'File', "T");
matlabFunction(T_s                 , 'Vars', {q, params, s}, 'File', "T_s");

%%
matlabFunction(v_i_1_i             , 'Vars', {q, dq, [L_0]}, 'File', "v_rel");
matlabFunction(omega_i_1_i         , 'Vars', {q, dq, [L_0]}, 'File', "omega_rel");
matlabFunction(dv_i_1_i            , 'Vars', {q, dq, ddq, [L_0]}, 'File', "a_rel");
matlabFunction(domega_i_1_i        , 'Vars', {q, dq, ddq, [L_0]}, 'File', "domega_rel");
matlabFunction(v_par_i             , 'Vars', {q, [L_0]}, 'File', "v_par");
matlabFunction(omega_par_i         , 'Vars', {q, [L_0]}, 'File', "omega_par");

cd(cF)
cd(destination_folder_b);
anonymousFunction(r_i_s, "./r_i_s");
anonymousFunction(dr_i_s, "./dr_i_s");
anonymousFunction(ddr_i_s, "./ddr_i_s");
anonymousFunction(J_dr_i_s, "./J_dr_i_s");
anonymousFunction(p_comi_s_a, "./p_comi_s");
anonymousFunction(dp_comi_s_a, "./dp_comi_s");
anonymousFunction(ddp_comi_s_a, "./ddp_comi_s");
%matlabFunction(p_comi              , 'Vars', {q, params}, 'File',"p_com");
anonymousFunction(p_comi, "./p_com");
%matlabFunction(v_i_comi            , 'Vars', {q, dq, params}, 'File', "v_com_rel");
anonymousFunction(v_i_comi, "./v_com_rel");
%matlabFunction(dv_i_comi           , 'Vars', {q, dq, ddq, params}, 'File', "a_com_rel");
anonymousFunction(dv_i_comi, "./a_com_rel");
anonymousFunction(I, "./I");
matlabFunction(m                   , 'Vars', {params}, 'File', "m");
matlabFunction(dI                  , 'Vars', {q, dq, params}, 'File', "dI");
anonymousFunction(J                , "./J");
matlabFunction(int_dr_i            , 'Vars', {q, dq, params}, 'File', "int_dr");
matlabFunction(int_ddr_i           , 'Vars', {q, dq, ddq,  params}, 'File', "int_ddr");
anonymousFunction(int_r_i_X_dr_i   , "./int_r_X_dr");
anonymousFunction(int_r_i_X_ddr_i  , "./int_r_X_ddr");
anonymousFunction(int_dr_i_X_pv_r  , "./int_dr_X_pv_r");
anonymousFunction(int_pv_r_O_dd_r  , "./int_pv_r_O_dd_r");
matlabFunction(int_dr_i_O_dr_i     , 'Vars', {q, dq,  params}, 'File', "int_dr_O_dr");matlabFunction(grad_int_dr_i       , 'Vars', {q,  params}, 'File', "grad_int_dr");
anonymousFunction(grad_int_r_i_X_dr_i , "./grad_int_r_X_dr");
anonymousFunction(grad_J_s         , "./grad_J_s");
anonymousFunction(grad_J           , "./grad_J");
anonymousFunction(J_r_i_X_dr_i_s   , "./J_r_i_X_dr_i_s");
%matlabFunction(grad_v_com_i        , 'Vars', {q,  params}, 'File', "grad_v_com");
anonymousFunction(grad_v_com_i     , "./grad_v_com");
matlabFunction(xi_i                , 'Vars', {q, params, s}, 'File', 'xi');
anonymousFunction(K_q              , "./K");
anonymousFunction(D_q              , "./D");
cd(cF)
