function logT = logmat(g)
%LOGMAT Logaritm map

theta = acos(trace(g)/2-1);

logT  = 1/8*csc(theta/2)^3*sec(theta/2)*((theta*cos(2*theta)-sin(theta))*eye(4) + ...
        - (theta*cos(theta) + 2*theta*cos(2*theta) - sin(theta) - sin(2*theta))*g + ...
        + (2*theta*cos(theta) + theta*cos(2*theta) - sin(theta) - sin(2*theta))*mpower(g, 2) + ...
        - (theta*cos(theta)-sin(theta))*mpower(g, 3));

end

