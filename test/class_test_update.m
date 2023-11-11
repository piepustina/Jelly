b1.Update(q_test, dq_test, ddq_test);
b2.updateBody(q_test, dq_test, ddq_test);
j2.updateJoint(q_test, dq_test, ddq_test);

%% Compute the transformation matrix of the tip
disp("Transformation matrix at tip...");
b1.T_
j2.T_

%% Compute the linear and angular velocity
disp("Linear tip velocity...");

b1.v_rel_
j2.v_rel_

disp("Body angular velocity...");
b1.omega_rel_
j2.omega_rel_

%% Compute the linear and angular acceleration
disp("Linear tip acceleration...");

b1.a_rel_
j2.a_rel_

disp("Body angular acceleration...");
b1.domega_rel_
j2.domega_rel_

%% Compute the CoM position
disp("CoM position...");
b1.p_com_
b2.p_com_

%% Compute the CoM velocity
disp("CoM velocity...");
b1.v_com_rel_
b2.v_com_rel_

%% Compute the CoM acceleration
disp("CoM acceleration...");
b1.a_com_rel_
b2.a_com_rel_

%% Compute the inertia
disp("Rotational inertia ...");
b1.I_
b2.I_

%% Compute the mass
disp("Body mass...");
b1.m_
b2.m_

%% Compute the time derivative of the inertia
disp("J...");
b1.J_
b2.J_

%% Compute int_dr
disp("int_dr...");
b1.int_dr_
b2.int_dr_

%% Compute int_ddr
disp("int_ddr...");
b1.int_ddr_
b2.int_ddr_

%% Compute int_r_X_dr
disp("int_r_X_dr...");
b1.int_r_X_dr_
b2.int_r_X_dr_

%% Compute int_r_X_ddr
disp("int_r_X_ddr...");
b1.int_r_X_ddr_
b2.int_r_X_ddr_

%% Compute int_dr_X_pv_r
disp("int_dr_X_pv_r...");
b1.int_dr_X_pv_r_
b2.int_dr_X_pv_r_

%% Compute int_dr_X_pv_r
disp("int_pv_r_O_dd_r ...");
b1.int_pv_r_O_dd_r_
b2.int_pv_r_O_dd_r_

%% Compute int_dr_X_pv_r
disp("grad_int_dr ...");
b1.grad_int_dr_
b2.grad_int_dr_

%% Compute grad_int_r_X_dr
disp("grad_int_r_X_dr ...");
b1.grad_int_r_X_dr_
b2.grad_int_r_X_dr_

%% Compute int_dr_O_dr
disp("int_dr_O_dr ...");
b1.int_dr_O_dr_
b2.int_dr_O_dr_

%% Compute int_dr_O_dr
disp("grad_int_dr ...");
b1.grad_int_dr_
b2.grad_int_dr_


%% Compute grad_J
disp("grad_J ...");
b1.grad_J_-b2.grad_J_

%% Compute grad_J
disp("grad_v_com ...");
b1.grad_v_com_
b2.grad_v_com_

%% Compute K
disp("K ...");
b1.K(q_test)
b2.K(q_test, Parameters)

%% Compute D
disp("D ...");
b1.D(q_test, dq_test)
b2.D(q_test, dq_test, Parameters)







