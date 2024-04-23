%#codegen
function s = sgn(x)
xSign       = x >= 0; 
s = double(xSign) -1.*double(~xSign);
end