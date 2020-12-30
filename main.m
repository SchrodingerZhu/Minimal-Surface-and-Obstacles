% dimensions
m = 6;
n = 6;
dim_x = (m-2) * (n-2);


% initial point
x0 = ones(dim_x, 1);

% boundary function
r = @(x,y) 0.5-abs(y-0.5);

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

% % C-BFGS
% c_opts = struct('m', 25, 'epsilon', 1e-5, 'eta', 0.2, 'gamma', 0.1, 'sigma', 0.5, 's', 1);
% L-BFGS
s_opts = struct('s', 1, 'sigma', 0.5, 'gamma', 0.1, 'm', 10, 'epsilon', 1e-7, 'delta', 1e-4, 'H', eye(dim_x)); % notice N is the size of X's vector


% minimization algorithm
% [x, obj, ~, ~] = gm_armijo(obj, grad, x0, opts);
% [ x, opt, G ] = L_BFGS (obj, grad, x0, s_opts);
[x,obj,iter] = momentum(obj,grad,x0,opts);