%% Build the body trees
r1 = BodyTree({j1}, {b1}, 1e-4*ones(n, 1), 1e-4*ones(n, 1), 1e-4*ones(n, 1), 1e-4*ones(n, 1));
r2 = BodyTree({j2}, {b2}, 1e-4*ones(n, 1), 1e-4*ones(n, 1), 1e-4*ones(n, 1), 1e-4*ones(n, 1));

disp("Mass matrix...")
r1.MassMatrix(q_test, 'double')
r2.MassMatrix(q_test, 'double')

disp("Stiffness force...")
r1.K(q_test, 'double')
r2.K(q_test, 'double')

disp("Gravity force...")
r1.GravityForce(q_test, 'double')
r2.GravityForce(q_test, 'double')

%% Convert the object to struct and back
sr1 = TreeStructConverter.ObjectToStruct(r1);
sr2 = TreeStructConverter.ObjectToStruct(r2);

r1n = TreeStructConverter.StructToObject(sr1);
r2n = TreeStructConverter.StructToObject(sr2);

r1n.MassMatrix(q_test, 'double')
r2n.MassMatrix(q_test, 'double')


