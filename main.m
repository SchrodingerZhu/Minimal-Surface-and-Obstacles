% dimensions
m = 6;
n = 6;
dim_x = (m-2) * (n-2);

% initial point
x0 = rand(dim_x, 1);

% boundary function
r = @(x,y) x^2+y^2;

% objective function and gradient
obj = @(x) f(x, m, n, r);
grad = @(x) g(x, m, n, r);

% options
opts.tol = 1e-5;
opts.s = 1;
opts.sigma = 0.5;
opts.gamma = 0.5;
opts.maxit = 100;
opts.maxit_armijo = 15;
opts.step = 1;
opts.alpha_lb = 1e-4;
opts.alpha_ub = 1e4;
opts.delta = 1;
opts.M = 3;

% minimization algorithm
[x, obj, ~, ~] = gm_bb_nonmonotone(obj, grad, x0, opts);