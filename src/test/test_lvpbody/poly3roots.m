clear; clc,
syms a x real
syms c real positive

p = (2*a/3)*x^3 + x^2 - c;

eqns = [p, a>0];

sol = solve(eqns, x, "ReturnConditions", true)