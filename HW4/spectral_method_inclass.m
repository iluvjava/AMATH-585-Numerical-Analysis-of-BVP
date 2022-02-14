d = domain(0, 1);

% Differential operator, linear or non-linear, with boundary conditions and
% the domain we are in. 
L = chebop(@(x, u) -diff(u, 2) + exp(x)*u, d, 0, 0);

% The RHS function. 
f = chebfun(@(x) - 3 + pi^2*sin(pi*x)+ exp(x)*(x^2 - x + sin(pi*x)), d);
u = L\f;

