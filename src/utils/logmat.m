function logT = logmat(g)
%LOGMAT Logaritm map

THETA_THSD = 1e-6;

ctheta = trace(g)/2-1;
% Make sure that the cos(theta) is in the interval [-1, 1]. In case map in
% into this interval.
if ctheta > 1
    ctheta = 1.0;
end
if ctheta < -1
    ctheta = -1.0;
end

theta = acos(ctheta);

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

