%% Build the body trees
r1 = BodyTree(J1, B1, 1e-4*ones(n*N_B, 1), 1e-4*ones(n*N_B, 1), 1e-4*ones(n*N_B, 1), 1e-4*ones(n*N_B, 1));
r2 = BodyTreeOld(J2, B2, 1e-4*ones(n*N_B, 1), 1e-4*ones(n*N_B, 1), 1e-4*ones(n*N_B, 1), 1e-4*ones(n*N_B, 1));

%%
disp("Mass matrix...")
r1.MassMatrix(repmat(q_test, N_B, 1), 'double')
r2.MassMatrix(repmat(q_test, N_B, 1), 'double')

disp("Stiffness force...")
r1.K(repmat(q_test, N_B, 1), 'double')
r2.K(repmat(q_test, N_B, 1), 'double')

disp("Gravity force...")
r1.GravityForce(repmat(q_test, N_B, 1), 'double')
r2.GravityForce(repmat(q_test, N_B, 1), 'double')