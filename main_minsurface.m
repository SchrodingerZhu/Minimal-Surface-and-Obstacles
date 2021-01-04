%% input data
% dimensions
m = 50;
n = 50;
dim_x = (m-2) * (n-2);

% initial point
rng(2021)
x0 = 2 * rand(dim_x,1);

% boundary function
% r = @(x,y) 1 + sin(2*pi*x);
% r = @(x,y) 1/2 - abs(y-1/2);
% r = @(x,y) 1 + asin(-1 + 2*sqrt(x*y));
% r = @(x,y) 1 + cos(1/(x+0.001));
% r = @(x,y) 1/(1+exp(x*y));
r = @(x,y) 5*x + (5*x) * sin(4*pi*x) -  y*cos(4*pi*y);
% r = @(x,y) -2*sin(2*pi*x) + 2*sin(2*pi*y);

% plot the initial surface
figure()
X0 = reshape(x0, [m-2, n-2]);
tri_visual (0, 1, 0, 1, addbd(X0, r));
title("Initial Surface")
% hold on

% objective function, gradient, and hessian
obj = @(x) objective(x, m, n, r);
grad = @(x) gradient(x, m, n, r);
hess = @(x) hessian(x, m, n, r);


%% options
% gradient method with Armijo backtracking line search
g_opts = struct("maxit", 1e5, "tol", 1e-5, "gamma", 0.5, "sigma", 0.5, "s", 1);

% gradient method with momentum
m_opts = struct("beta", 0.5, "tol", 1e-5, "L0", 0.1, "sigma", 0.5);

% gradient method with Barzilai-Borwein step
bb_opts = struct("maxit", 1e4, "tol", 1e-5, "step", 2);

% gradient method with Barzilai-Borwein step and nonmonotone line search
bbn_opts = struct("maxit", 1e4, "tol", 1e-5, "alpha_lb", 1e-4, "alpha_ub", 1e4, "delta", 1, "gamma", 0.5, "sigma", 0.5, "M", 5);

% globalized Newton"s method
newton_opts = struct("maxit", 1e3, "tol", 1e-5, "beta1", 1e-6, "beta2", 1e-6, "p", 0.1, "gamma", 0.4, "sigma", 0.5, "s", 1);

% L-BFGS
s_opts = struct("s", 1, "sigma", 0.5, "gamma", 0.1, "m", 25, "epsilon", 1e-5, "eta", 0.5, "H", eye(dim_x)); 

% C-BFGS
c_opts = struct("m", 25, "epsilon", 1e-5, "eta", 0.2, "gamma", 0.1, "sigma", 0.5, "s", 1);


%% minimization algorithms
% [x, opt_val, f_k, grad_k, T, method] = gm_armijo(obj, grad, x0, g_opts);
% [x,opt_val, f_k, grad_k, T, method] = gm_momentum(obj, grad, x0, m_opts);
% [x, opt_val, f_k, grad_k, T, method] = gm_bb(obj, grad, x0, bb_opts);
% [x, opt_val, f_k, grad_k, T, method] = gm_bb_nonmonotone(obj, grad, x0, bbn_opts);
[x, opt_val, f_k, grad_k, T, method] = newton_globalized(obj, grad, hess, x0, newton_opts);
% [x, opt_val, f_k, grad_k, T, method] = PW_L_BFGS (obj, grad, x0, s_opts);
% [x, opt_val, f_k, grad_k, T, method] = C_BFGS (obj, grad, x0, c_opts);
fprintf("CPU Time = %.4f\n", T(end))

%% Plots
% plot the generated minimal surface
figure()
X = reshape(x, [m-2, n-2]);
tri_visual (0, 1, 0, 1, addbd(X, r));
title("Minimal Surface", "interpreter", "latex")

% plot iterates
plot_iter(f_k, grad_k, T, method)