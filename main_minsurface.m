%% input data
% dimensions
m = 12;
n = 12;
dim_x = (m-2) * (n-2);

% initial point
x0 = 2 * rand(dim_x,1);

% boundary function
r = @(x,y) 1 + sin(2*pi*x);

% plot the initial surface
figure(1)
X0 = reshape(x0, [m-2, n-2]);
tri_visual (0, 1, 0, 1, addbd(X0, r)); 
% hold on

% objective function, gradient, and hessian
obj = @(x) objective(x, m, n, r);
grad = @(x) gradient(x, m, n, r);
hess = @(x) hessian(x, m, n, r);

%% symbolic computation of gradient and hessian
% y = sym('y', [(m-2)*(n-2),1]);
% F = obj(y);
% 
% % gradient
% fprintf("===== Gradient =====\n")
% fprintf("- - - differentiation begins\n")
% grad = gradient(F, y);
% hess = jacobian(grad, y);
% fprintf("- - - differentiation ends\n")
% 
% fprintf("- - - coverting begins\n")
% % grad = @(x) g(x, m, n, r);
% grad = matlabFunction(grad, "vars", {y});
% hess = matlabFunction(hess, "vars", {y});
% fprintf("- - - converting ends\n")

%% options
% gradient method with Armijo backtracking line search
g_opts = struct('maxit', 1e4, 'tol', 1e-5, 'gamma', 0.5, 'sigma', 0.5, 's', 1);

% gradient method with Barzilai-Borwein step
bb_opts = struct('maxit', 1e3, 'tol', 1e-5, 'step', 1);

% gradient method with Barzilai-Borwein step and nonmonotone line search
bbn_opts = struct('maxit', 1e3, 'tol', 1e-5, 'alpha_lb', 1e-4, 'alpha_ub', 1e4, 'delta', 1, 'gamma', 0.5, 'sigma', 0.5, 'M', 5);

% globalized Newton's method
newton_opts = struct('maxit', 1e3, 'tol', 1e-5, 'beta1', 1e-6, 'beta2', 1e-6, 'p', 0.1, 'gamma', 0.4, 'sigma', 0.5, 's', 1);

% C-BFGS
c_opts = struct('m', 25, 'epsilon', 1e-5, 'eta', 0.2, 'gamma', 0.1, 'sigma', 0.5, 's', 1);

% L-BFGS
s_opts = struct('s', 1, 'sigma', 0.5, 'gamma', 0.1, 'm', 10, 'epsilon', 1e-7, 'delta', 1e-4, 'H', eye(dim_x)); % notice N is the size of X's vector


%% minimization algorithms
% [x, opt_val, ~, ~] = gm_armijo(obj, grad, x0, g_opts);
% [x, opt_val, ~, ~] = gm_bb(obj, grad, x0, bb_opts);
% [x, opt_val, ~, ~] = gm_bb_nonmonotone(obj, grad, x0, bbn_opts);
[x,opt_val,f_k,grad_k,h_paul] = newton_globalized(obj, grad, hess, x0, newton_opts);
% [x, opt_val, G] = L_BFGS (obj, grad, x0, s_opts);
% [x, opt_val, G] = C_BFGS (obj, grad, x0, c_opts);
% [x,opt_val,iter] = momentum(obj,grad,x0,opts);
toc

% plot the generated minimal surface
figure("After Optimization")
X = reshape(x, [m-2, n-2]);
tri_visual (0, 1, 0, 1, addbd(X, r)); 