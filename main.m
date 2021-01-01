%% input data
% dimensions
m = 30;
n = 30;
dim_x = (m-2) * (n-2);

% initial point
x0 = rand(dim_x, 1);

% boundary function
r = @(x,y) 1+cos(1/(x+0.001));

figure(1)
X0 = reshape(x0, [m-2, n-2]);
tri_visual (0, 1, 0, 1, addbd(X0, r)); 
% hold on
% objective function and gradient
obj = @(x) f(x, m, n, r);
grad = @(x) g(x, m, n, r);

%% options
% gradient method with Armijo backtracking line search
g_opts = struct('maxit', 1e4, 'tol', 1e-5, 'gamma', 0.5, 'sigma', 0.5, 's', 1);

% gradient method with Barzilai-Borwein step
bb_opts = struct('maxit', 1e3, 'tol', 1e-5, 'step', 1);

% gradient method with Barzilai-Borwein step and nonmonotone line search
bbn_opts = struct('maxit', 1e3, 'tol', 1e-5, 'alpha_lb', 1e-4, 'alpha_ub', 1e4, 'delta', 1, 'gamma', 0.5, 'sigma', 0.5, 'M', 5);

% C-BFGS
c_opts = struct('m', 25, 'epsilon', 1e-5, 'eta', 0.2, 'gamma', 0.1, 'sigma', 0.5, 's', 1);
% L-BFGS
s_opts = struct('s', 1, 'sigma', 0.5, 'gamma', 0.1, 'm', 10, 'epsilon', 1e-7, 'delta', 1e-4, 'H', eye(dim_x)); % notice N is the size of X's vector


%% minimization algorithms
% [x, obj, ~, ~] = gm_armijo(obj, grad, x0, g_opts);
[x, obj, ~, ~] = gm_bb(obj, grad, x0, bb_opts);
% [x, obj, ~, ~] = gm_bb_nonmonotone(obj, grad, x0, bbn_opts);
% [x, opt, G] = L_BFGS (obj, grad, x0, s_opts);
% [x, opt, G] = C_BFGS (obj, grad, x0, c_opts);
% [x,obj,iter] = momentum(obj,grad,x0,opts);

% plot
figure(2)
X = reshape(x, [m-2, n-2]);
tri_visual (0, 1, 0, 1, addbd(X, r)); 