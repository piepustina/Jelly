%% Compute the transformation matrix of the tip
disp("Transformation matrix at tip...");
b1.T(q_test)
j1.T(q_test)
j2.T(q_test, Parameters(1))


%% Compute transformation matrix at half rest length
disp("Transformation matrix at s = L0/2");
b1.T_s(q_test, L0/2)
j2.T_s(q_test, Parameters(1), L0/2)

%% Compute the linear and angular velocity
disp("Linear tip velocity...");

b1.v_rel(q_test, dq_test)
j1.v_rel(q_test, dq_test)
j2.v_rel(q_test, dq_test, Parameters(1))

disp("Body angular velocity...");
b1.omega_rel(q_test, dq_test)
j1.omega_rel(q_test, dq_test)
j2.omega_rel(q_test, dq_test, Parameters(1))

%% Compute the linear and angular acceleration
disp("Linear tip acceleration...");

b1.a_rel(q_test, dq_test, ddq_test)
j1.a_rel(q_test, dq_test, ddq_test)
j2.a_rel(q_test, dq_test, ddq_test, Parameters(1))

disp("Body angular acceleration...");
b1.domega_rel(q_test, dq_test, ddq_test)
j1.domega_rel(q_test, dq_test, ddq_test)
j2.domega_rel(q_test, dq_test, ddq_test, Parameters(1))

%% Compute the linear and angular Jacobians
disp("Linear Jacobian...");

b1.v_par(q_test)
j1.v_par(q_test)
j2.v_par(q_test, Parameters(1))

disp("Angular jacobian...");
b1.omega_par(q_test)
j1.omega_par(q_test)
j2.omega_par(q_test, Parameters(1))

%% Compute the CoM position
disp("CoM position...");
b1.p_com(q_test)
b2.p_com(q_test, Parameters)

%% Compute the CoM velocity
disp("CoM velocity...");
b1.v_com_rel(q_test, dq_test)
b2.v_com_rel(q_test, dq_test, Parameters)

%% Compute the CoM acceleration
disp("CoM acceleration...");
b1.a_com_rel(q_test, dq_test, ddq_test)
b2.a_com_rel(q_test, dq_test, ddq_test, Parameters)

%% Compute the inertia
disp("Rotational inertia ...");
b1.I(q_test)
b2.I(q_test, Parameters)

%% Compute the mass
disp("Body mass...");
b1.m()
b2.m(Parameters)

%% Compute the time derivative of the inertia
disp("J...");
b1.J(q_test, dq_test)
b2.J(q_test, dq_test, Parameters)

%% Compute int_r_X_dr
disp("int_r_X_dr...");
b1.int_r_X_dr(q_test, dq_test)
b2.int_r_X_dr(q_test, dq_test, Parameters)

%% Compute int_r_X_ddr
disp("int_r_X_ddr...");
b1.int_r_X_ddr(q_test, dq_test, ddq_test)
b2.int_r_X_ddr(q_test, dq_test, ddq_test, Parameters)

%% Compute int_dr_X_pv_r
disp("int_dr_X_pv_r...");
b1.int_dr_X_pv_r(q_test, dq_test)
b2.int_dr_X_pv_r(q_test, dq_test, Parameters)

%% Compute int_dr_X_pv_r
disp("int_pv_r_O_dd_r ...");
b1.int_pv_r_O_dd_r(q_test, dq_test, ddq_test)
b2.int_pv_r_O_dd_r(q_test, dq_test, ddq_test, Parameters)

%% Compute int_dr_X_pv_r
disp("grad_int_dr ...");
b1.grad_int_dr(q_test)
b2.grad_int_dr(q_test, Parameters)

%% Compute grad_int_r_X_dr
disp("grad_int_r_X_dr ...");
b1.grad_int_r_X_dr(q_test)
b2.grad_int_r_X_dr(q_test, Parameters)


%% Compute grad_J
disp("grad_J ...");
b1.grad_J(q_test)
b2.grad_J(q_test, Parameters)

%% Compute grad_J
disp("grad_v_com ...");
b1.grad_v_com(q_test)
b2.grad_v_com(q_test, Parameters)

%% Compute K
disp("K ...");
b1.K(q_test)
b2.K(q_test, Parameters)

%% Compute D
disp("D ...");
b1.D(q_test, dq_test)
b2.D(q_test, dq_test, Parameters)







