% Test the plot functionality of the soft robot class
%Create an instance of the class
clear; clc;
L0              = 0.5;
BaseRadius      = 0.01;
TipRadius       = 0.001;
MassDensity     = 1062;
YoungModulus    = 6.66e5;
PoissonRatio    = 0.5;
MaterialDamping = 0.1;
NGaussPoints    = 10;
N_B             = 1;
Parameters      = [L0, BaseRadius, TipRadius, MassDensity, YoungModulus, PoissonRatio, MaterialDamping]';
B1              = cell(N_B, 1);
J1              = cell(N_B, 1);
for i = 1:N_B
    B1{i} = PCC3D([Parameters; NGaussPoints]);
    J1{i} = FixedJoint();
end

%Build the soft robot
sr = SoftRobot(J1, B1);

%%
figure; hold on; grid on; 
light; lighting gouraud
%Generate a random configuration
q_lb   = repmat([-pi; -pi; -0.05], N_B, 1);
q_ub   = repmat([ pi;  pi;  0.05], N_B, 1);
q_test = q_lb + (q_ub-q_lb).*rand(sr.n,1);
sr.plot(q_test);

view(3)
axis equal



