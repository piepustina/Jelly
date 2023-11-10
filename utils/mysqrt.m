%#codegen
function sqrtX = mysqrt(X)
%MYSQRT Implements the "square" root for X real but possibily negative.
%This is supposed to be used only in combination with the fresnel integrals
% if X < 0
%     sqrtX = sqrt(-X);
% else
%     sqrtX = sqrt(X);
% end
sqrtX = complex(0, 0);
if X < 0
    sqrtX = sqrt(complex(X));
else
    sqrtX = sqrt(X);
end
end

