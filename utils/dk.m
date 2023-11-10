function tip_pos = dk(q, L, A0)

n = length(q);

T_prev = [A0(1:3,1:2) A0(2:4,4)];

for i = 1 : n
    %Check if the curvature is infinite
    if abs(q(i)) < 0.0017%set a threshold for the infinite curvature (below 0.1 deg)
        T = T_prev*[cos(q(i)) -sin(q(i)) L(i);
                    sin(q(i))  cos(q(i)) 0;
                    0          0         1];
    else
        %Compute the transformation matrix from the current frame to the base
        T = T_prev*[cos(q(i)) -sin(q(i)) L(i)*sin(q(i))/q(i);
                    sin(q(i))  cos(q(i)) L(i)*(1 - cos(q(i)))/q(i);
                    0          0         1];
    end
    %Prepare for the next iteration
    T_prev = T;
end

tip_pos = T(1:2, 3);


end

