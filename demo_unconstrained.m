% objective function and gradient
f = @(x) (exp(1-x(1)-x(2)) + exp(x(1)+x(2)-1) + x(1)^2 + x(1)*x(2) + x(2)^2 + 2*x(1) - 3*x(2));
fm = @(x1, x2) (exp(1-x1-x2) + exp(x1+x2-1) + x1^2 + x1*x2 + x2^2 + 2*x1 - 3*x2);   % multivatiable form
grad1 = @(x) (-exp(1-x(1)-x(2)) + exp(x(1)+x(2)-1) + 2*x(1) + x(2) + 2);
grad2 = @(x) (-exp(1-x(1)-x(2)) + exp(x(1)+x(2)-1) + 2*x(2) + x(1) - 3);
grad = @(x) [grad1(x); grad2(x)];

% parameters
x0 = [0; 0];   % initial point
opts.tol = 1e-5;   % tolerance
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

[x,obj,f_k,grad_k] = gm_armijo(f, grad, x0, opts);
for i = 1 : 2
    opts.step = i;
    [x,obj,f_k,grad_k] = gm_bb(f, grad, x0, opts);
end
[x,obj,f_k,grad_k] = gm_bb_nonmonotone(f, grad, x0, opts);