function Aq = symbolicActuationMatrix(Tree, q, Start, End)
%NOTE: The method assumes that the actuators start at the base of a body
%and terminate at the end of another.

%SYMBOLICACTUATIONMATRIX Given a BodyTree Tree computes the actuation
%matrix is symbolic form, using q as configuration variables and s as abscissa curvlinear.
%This method works also when the Tree contains
%symbolic terms. Start is a vector such that its i-th element indicate the
%Body where the actuator starts/ends.

%Curvilinear abscissa
s = sym('s', [1, 1], 'real');

%Symbolic expression for the actuation matrix
Aq = zeros(Tree.n, Tree.N_A, 'sym');

for i = 1:Tree.N_A
    %Get the current actuator and its interval
    Act = Tree.Actuators{i};
    L_1 = Act.Interval(1);
    %Location of the actuator as a function of the curvilinear abscissa s
    d = Act.d(s);
    %Derivative w.r.t. s
    d_prime = diff(d, s);

    idx = 1;
    for j = 1:Tree.N_B
        if (j >= Start(i)) && (j <= End(i))
            Bj  = Tree.Bodies{j};
            q_j = q(idx:idx+Bj.n-1);
            xi = Bj.xi(q_j, Bj.Parameters, s);

            %Length of the current body
            L_2 = L_1 + Bj.RestLength;
            
            %Compute the unit tangent vector to each actuator
            xi_skew = skew(xi(1:3));
            tc = xi_skew*d + xi(4:6) + d_prime;
            t  = (tc./sqrt(tc'*tc));
            
            %Compute the column of the actuation matrix
            d_skew = skew(d);
            Phi_a = [d_skew*t; t];
            xi_par = jacobian(xi, q_j);
            Aq(idx:idx+Bj.n-1, i) = int(xi_par'*Phi_a, s, [L_1, L_2]);
    
            %Prepare for the next iteration
            L_1 = L_2;
        end

        idx = idx+Tree.Bodies{j}.n;
    end
            
end

end

