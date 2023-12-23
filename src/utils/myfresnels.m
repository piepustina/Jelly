function fX = myfresnels(X)
%MYFRESNELS Computes the Fresnel sine at X with support to code generation
% The function is implemented using eq.(13) of https://www.hindawi.com/journals/mpe/2018/4031793/
wasComplex = 0;
if ~isreal(X)
    X = imag(X);
    wasComplex = 1;
end
fX = sign(X)*(1/2 - (1+0.926*X)/(2+1.792*X+3.104*X^2)*cos(pi*X^2/2)-...
           1/(2+4.142*X+3.492*X^2+6.67*X^3)*sin(pi*X^2/2));
if wasComplex
    fX = -1i*fX;
end
end