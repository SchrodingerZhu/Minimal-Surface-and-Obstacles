% demo
beta = 0.75;
A = [1 0 1; 0 1 0; 1 0 1];
f = @(x) x'*A*x;
h = @(x) x'*x - 1;

grad_f = @(x) 2*A*x;
grad_h = @(x) 2*x;

eps = @(k) 1e-3 / k;
alph = @(k) 1;

x0 = [1;0;0];
mu0 = 0.5;

opts.maxit = 100;
opts.tol = eps;
opts.alpha = alph;

[x, obj] = augment_lagrange(f, h, grad_f, grad_h, x0, mu0, 2, opts);