function lvp_codegen(Nodes, Elements)
%LVP_CODEGEN Test code generation for the locally volume preserving body

coder.gpu.kernelfun;

Primitives = cell(1, LVPBody.MaxPrimitivesNumber);
Primitives{1} = PCStretchCompressionPrimitive(0.3);
Primitives{2} = PCTwistShearPrimitive(0.3);
Primitives{3} = PCBendingPrimitive(0.3);
for i = 4:LVPBody.MaxPrimitivesNumber
    Primitives{i} = 0;
end

tic, B = LVPBody(Nodes, Elements, Primitives, 2); toc;

q_test = [2; 1; -1; 1; 2; 1];

tic, B.UpdateKinematics(q_test, 0*q_test, 0*q_test); toc;

tic, B.UpdateKinematics(2*q_test, 0*q_test, 0*q_test); toc;

tic, B.Update(q_test, 0*q_test, 0*q_test); toc;

tic, B.Update(2*q_test, 0*q_test, 0*q_test); toc;

end

