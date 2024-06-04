function logT = logmat(T)
%LOGMAT Logaritm map of T, where T is an homogeneous transformation matrix or a 3D tensor of transformation matrices.

arguments (Input)
    T (4, 4, :)
end

arguments (Output)
    logT (4, 4, :)
end

% Threshold for the numerical evaluation of the logarithm map.
THETA_THSD = 1e-6;

nT = size(T, 3);

% Cosine of theta
ctheta = (T(1, 1, 1:nT) + T(2, 2, 1:nT) + T(3, 3, 1:nT))/2-1;
% Make sure that the cos(theta) is in the interval [-1, 1]. In case map in
% into this interval.
ctheta(ctheta > 1)  = 1.0;
ctheta(ctheta < -1) = -1.0;
theta               = acos(ctheta);

% Approximate theta around zero with the threshold value
ThetaIdx            = abs(theta) <= THETA_THSD;
theta(ThetaIdx)     = sign(theta(ThetaIdx)).*THETA_THSD;
theta(theta == 0)   = THETA_THSD;

% Compute the output
ctheta  = cos(theta);
stheta  = sin(theta);
c2theta = cos(2*theta);
s2theta = sin(2*theta);

powT2   = pagemtimes(T, T);
powT3   = pagemtimes(powT2, T);

% Non vectorized form
% logT  = 1/8*csc(theta/2)^3*sec(theta/2)*((theta*cos(2*theta)-sin(theta))*eye(4) + ...
%         - (theta*cos(theta) + 2*theta*cos(2*theta) - sin(theta) - sin(2*theta))*T + ...
%         + (2*theta*cos(theta) + theta*cos(2*theta) - sin(theta) - sin(2*theta))*mpower(T, 2) + ...
%         - (theta*cos(theta)-sin(theta))*mpower(T, 3));

logT  = pagemtimes(1/8*(csc(theta/2).^3).*sec(theta/2), ...
        pagemtimes((theta.*c2theta-stheta), repmat(eye(4), 1, 1, nT)) ...
      - pagemtimes(theta.*ctheta + 2*theta.*c2theta - stheta - s2theta, T)...
      + pagemtimes(2*theta.*ctheta + theta.*c2theta - stheta - s2theta, powT2)...
      - pagemtimes(theta.*ctheta-stheta, powT3));
end

