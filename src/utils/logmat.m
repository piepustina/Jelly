function logT = logmat(g)
%LOGMAT Logaritm map

THETA_THSD = 1e-6;

theta = acos(trace(g)/2-1);

if abs(theta) <= THETA_THSD
    theta = sign(theta)*THETA_THSD;
end

if theta == 0
    theta = THETA_THSD;
end

logT  = 1/8*csc(theta/2)^3*sec(theta/2)*((theta*cos(2*theta)-sin(theta))*eye(4) + ...
        - (theta*cos(theta) + 2*theta*cos(2*theta) - sin(theta) - sin(2*theta))*g + ...
        + (2*theta*cos(theta) + theta*cos(2*theta) - sin(theta) - sin(2*theta))*mpower(g, 2) + ...
        - (theta*cos(theta)-sin(theta))*mpower(g, 3));

end

