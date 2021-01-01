% instance 1
% f = @(x) (exp(1-x(1)-x(2)) + exp(x(1)+x(2)-1) + x(1)^2 + x(1)*x(2) + x(2)^2 + 2*x(1) - 3*x(2));
% fm = @(x1, x2) (exp(1-x1-x2) + exp(x1+x2-1) + x1^2 + x1*x2 + x2^2 + 2*x1 - 3*x2);   % multivatiable form
% grad1 = @(x) (-exp(1-x(1)-x(2)) + exp(x(1)+x(2)-1) + 2*x(1) + x(2) + 2);
% grad2 = @(x) (-exp(1-x(1)-x(2)) + exp(x(1)+x(2)-1) + 2*x(2) + x(1) - 3);
% grad = @(x) [grad1(x); grad2(x)];

% instance 2
f = @(x) 100 * (x(2) - x(1)^2)^2 + (x(1)-1)^2;
grad = @(x) [400*(x(1)^2-x(2))*x(1) + 2*(x(1)-1); 200*(x(2)-x(1)^2)];
hess = @(x) [1200*x(1)^2-400*x(2)+2 -400*x(1); -400*x(1) 200];

% initial point
x0 = [0; 0];

% gradient method with Armijo backtracking line search
g_opts = struct('maxit', 1e4, 'tol', 1e-5, 'gamma', 0.5, 'sigma', 0.5, 's', 1);

% gradient method with Barzilai-Borwein step
bb_opts = struct('maxit', 1e3, 'tol', 1e-5, 'step', 1);

% gradient method with Barzilai-Borwein step and nonmonotone line search
bbn_opts = struct('maxit', 1e3, 'tol', 1e-5, 'alpha_lb', 1e-4, 'alpha_ub', 1e4, 'delta', 1, 'gamma', 0.5, 'sigma', 0.5, 'M', 5);

% globalized Newton's method
newton_opts = struct('maxit', 1e3, 'tol', 1e-5, 'beta1', 1e-6, 'beta2', 1e-6, 'p', 0.1, 'gamma', 0.4, 'sigma', 0.5, 's', 1);

% optimization
[x,obj,f_k,grad_k] = gm_armijo(f, grad, x0, g_opts);
for i = 1 : 2
    bb_opts.step = i;
    [x,obj,f_k,grad_k] = gm_bb(f, grad, x0, bb_opts);
end
[x,obj,f_k,grad_k] = gm_bb_nonmonotone(f, grad, x0, bbn_opts);
[x,obj,f_k,grad_k] = newton_globalized(f, grad, hess, x0, newton_opts);