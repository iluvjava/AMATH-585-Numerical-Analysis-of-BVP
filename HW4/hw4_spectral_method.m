d = domain(0, 1); 
L = chebop(@(x, u) -diff((1 + x^2)*diff(u, 1), 1), d, 0, 0); 
f = chebfun(@(x) 2*(3*x^2 - x + 1), d);
u = L\f;
plot(u, '-o');
x = linspace(0, 1, 10);
hold on;
plot(x, x.*(1 - x), 'x');
title("Cheb Spectral Solution");
xlabel("x");
ylabel("u(x)");
legend(["u(x)", "x(1 - x)"]);

correctU = @(x) x.*(1 - x);
infError = norm(correctU(x) - u(x), inf);
l2Error = sqrt(1/2)*norm(correctU(x) - u(x));
disp(infError);
disp(l2Error);