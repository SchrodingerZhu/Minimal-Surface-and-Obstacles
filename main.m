% dimensions
m = 100;
n = 100;
dim_x = (m-2) * (n-2);
% rng(2020);

% initial point
% x0 = [2.5*rand(dim_x/4,1) ; 2.2*rand(dim_x/4,1) ; 2*rand(dim_x/2,1)];
% x0 = x0(randperm(dim_x));
% x0 = ones(dim_x,1);
x0 = 20*rand(dim_x,1);

% boundary function
r = @(x,y) x+10;
% r = @(x,y) 0.5-abs(0.5-y);

figure(1)
X0 = reshape(x0, [m-2, n-2]);
[~] = tri_visual (0, 1, 0, 1, addbd(X0, r)); 

% objective function and gradient
obj = @(x) f(x, m, n, r);
grad = @(x) g(x, m, n, r);

% options
opts.tol = 1e-5;
opts.s = 1;
opts.sigma = 0.5;
opts.gamma = 0.5;
opts.maxit = 500;
opts.maxit_armijo = 15;
opts.step = 1;
opts.alpha_lb = 1e-4;
opts.alpha_ub = 1e4;
opts.delta = 1;
opts.M = 3;

% % C-BFGS
c_opts = struct('m', 25, 'epsilon', 1e-5, 'eta', 0.2, 'gamma', 0.1, 'sigma', 0.5, 's', 1);
% L-BFGS
s_opts = struct('s', 1, 'sigma', 0.5, 'gamma', 0.1, 'm', 10, 'epsilon', 1e-5, 'delta', 1e-4, 'H', eye(dim_x)); % notice N is the size of X's vector


% minimization algorithm
fprintf("before: f(x) = %.4f\n", obj(x0))
[x, obj_val, ~, ~] = gm_bb(obj, grad, x0, opts);
fprintf("after: f(x) = %.4f\n", obj(x))
% [ x, opt, G ] = C_BFGS (obj, grad, x0, c_opts);
% [x,obj,iter] = momentum(obj,grad,x0,opts);

% plot
figure(2)
X = reshape(x, [m-2, n-2]);
[~] = tri_visual (0, 1, 0, 1, addbd(X, r)); 