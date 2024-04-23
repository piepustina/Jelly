%#codegen
function Ts = Tinterpolate(T1, T2, s)
%TINTERPOLATE Compute the linear interpolation between two transformation matrices at the normalized point s
% s varies in the interval [0, 1]

%TODO: It is possible to speed up this function by passing 3D tensors for T1, T2 and s.

arguments (Input)
    T1 (4, 4) double
    T2 (4, 4) double
    s  (1, :) double
end

arguments (Output)
    Ts (4, 4, :) double
end

% Output preallocation
nS = length(s);
Ts = repmat(eye(4), 1, 1, nS);

% Linearly interpolate the position
Ts(1:3, 4, 1:nS) = T1(1:3, 4).*(1-s) + T2(1:3, 4).*s;


% Interpolate the rotational part using the slerp method
Quat1 = RotmToQuat(T1(1:3, 1:3));% Convert the start rotation into a quaternion
Quat2 = RotmToQuat(T2(1:3, 1:3));% Convert the final rotation into a quaternion

% Compute the dot product bewteen the quaternions
d                           = dot(Quat1, Quat2, 1);
if d < 0
    dotQuat1                = -Quat1;
else
    dotQuat1                =  Quat1;
end
Theta                       = acos(abs(d));

if abs(Theta) < 1e-4% Use the limit expression close to zero
    QuatInterp = (dotQuat1.*(1-s) + Quat2.*s);
    QuatInterp = QuatInterp./vecnorm(QuatInterp, 2);
else
    QuatInterp = (sin((1-s).*Theta).*dotQuat1 + sin(s.*Theta).*Quat2)./sin(Theta);
end

% Convert back to rotation
Rs = QuatToRotm(QuatInterp);
% Store the result into the output
Ts(1:3, 1:3, 1:nS) = Rs(1:3, 1:3, 1:nS);

if ~isreal(Ts)
    error("Errore non reale!");
end

end

