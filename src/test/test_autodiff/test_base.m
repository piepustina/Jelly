% Create a dual variable

x = Dual([2; 0]);

y = Dual([0; 0]);


x + y

%%
f = @(t) t^2;

f_res = f(x);

g = @(t) f_res*t;

