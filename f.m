function y = f(X, m, n) %make the function calculation efficient
y = sum(vecnorm(manip(addbd(X, m, n), m, n)));






