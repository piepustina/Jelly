function Expr = replaceFresnel(expr)
%REPLACEFRESNEL Replace the fresnel integrals from the input expression

disp("Replacing Fresnels");

Expr_string = string(expr);
Expr_string = strrep(Expr_string,'fresnelc','myfresnelc');
Expr_string = strrep(Expr_string,'fresnels','myfresnels');
Expr     = str2sym(Expr_string);

end

