%% Export all the functions for the current model
%--------------------------------------------------------------------------
% Store all the terms describing the body also in a common folder ---------
%--------------------------------------------------------------------------
disp("Exporting functions...");

if ~exist('optimize_export','var')
    optimize_export = true;
end
cF = pwd;
cd(destination_folder_j);

disp("Saving the Joint functions...");

matlabFunction(T                   , 'Vars', {q, [L_0]}, 'File', "T", 'Optimize', optimize_export);
matlabFunction(T_s                 , 'Vars', {q, params, s}, 'File', "T_s", 'Optimize', optimize_export);
matlabFunction(v_i_1_i             , 'Vars', {q, dq, [L_0]}, 'File', "v_rel", 'Optimize', optimize_export);
matlabFunction(omega_i_1_i         , 'Vars', {q, dq, [L_0]}, 'File', "omega_rel", 'Optimize', optimize_export);
matlabFunction(dv_i_1_i            , 'Vars', {q, dq, ddq, [L_0]}, 'File', "a_rel", 'Optimize', optimize_export);
matlabFunction(domega_i_1_i        , 'Vars', {q, dq, ddq, [L_0]}, 'File', "domega_rel", 'Optimize', optimize_export);
matlabFunction(v_par_i             , 'Vars', {q, [L_0]}, 'File', "v_par", 'Optimize', optimize_export);
matlabFunction(omega_par_i         , 'Vars', {q, [L_0]}, 'File', "omega_par", 'Optimize', optimize_export);

cd(cF)
cd(destination_folder_b);

disp("Saving the Body functions...");

matlabFunction(p_comi              , 'Vars', {q, params}, 'File',"p_com", 'Optimize', optimize_export);
matlabFunction(v_i_comi            , 'Vars', {q, dq, params}, 'File', "v_com_rel", 'Optimize', optimize_export);
matlabFunction(dv_i_comi           , 'Vars', {q, dq, ddq, params}, 'File', "a_com_rel", 'Optimize', optimize_export);
matlabFunction(I                   , 'Vars', {q, params}, 'File', "I", 'Optimize', optimize_export);
matlabFunction(m                   , 'Vars', {params}, 'File', "m", 'Optimize', optimize_export);
matlabFunction(dI                  , 'Vars', {q, dq, params}, 'File', "dI", 'Optimize', optimize_export);
matlabFunction(J                   , 'Vars', {q, dq, params}, 'File', "J", 'Optimize', optimize_export);
matlabFunction(int_dr_i            , 'Vars', {q, dq, params}, 'File', "int_dr", 'Optimize', optimize_export);
matlabFunction(int_ddr_i           , 'Vars', {q, dq, ddq,  params}, 'File', "int_ddr", 'Optimize', optimize_export);
matlabFunction(int_r_i_X_dr_i      , 'Vars', {q, dq,  params}, 'File', "int_r_X_dr", 'Optimize', optimize_export);
matlabFunction(int_r_i_X_ddr_i     , 'Vars', {q, dq, ddq,  params}, 'File', "int_r_X_ddr", 'Optimize', optimize_export);
matlabFunction(int_dr_i_X_pv_r     , 'Vars', {q, dq,  params}, 'File', "int_dr_X_pv_r", 'Optimize', optimize_export);
matlabFunction(int_pv_r_O_dd_r     , 'Vars', {q, dq, ddq,  params}, 'File', "int_pv_r_O_dd_r", 'Optimize', optimize_export);
matlabFunction(int_dr_i_O_dr_i     , 'Vars', {q, dq,  params}, 'File', "int_dr_O_dr", 'Optimize', optimize_export);
matlabFunction(grad_int_dr_i       , 'Vars', {q,  params}, 'File', "grad_int_dr", 'Optimize', optimize_export);
matlabFunction(grad_int_r_i_X_dr_i , 'Vars', {q,  params}, 'File', "grad_int_r_X_dr", 'Optimize', optimize_export);
matlabFunction(grad_J              , 'Vars', {q,  params}, 'File', "grad_J", 'Optimize', optimize_export);
matlabFunction(grad_v_com_i        , 'Vars', {q,  params}, 'File', "grad_v_com", 'Optimize', optimize_export);
matlabFunction(xi_i                , 'Vars', {q, params, s}, 'File', 'xi', 'Optimize', optimize_export);
matlabFunction(K_q                 , 'Vars', {q, params}, 'File', 'K', 'Optimize', optimize_export);
matlabFunction(D_q                 , 'Vars', {q, dq, params}, 'File', 'D', 'Optimize', optimize_export);
cd(cF)
