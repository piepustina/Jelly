function [X,Y] = buildSegment(q, L, n)
%buildSegment Function that computes a set of points [X,Y] representing a constant
%curvature segment of lenght L and with degree of curvature q. The segment
%is discretized with n points

dqs = linspace(0, q, n);
%If the curvature is infinite plot just a line
if abs(q) < 0.0017%set a threshold for the infinite curvature (below 0.1 deg)
    X = linspace(0, L, n);
    Y = linspace(0, 0, n);
else
    %The radius of curvature is not infinite, compute it.
    rho = L/q;
    %Compute the position of n equally spaced points along the CC segment
    X = rho.* sin(dqs);
    Y = rho.* (1-cos(dqs));
end

end

