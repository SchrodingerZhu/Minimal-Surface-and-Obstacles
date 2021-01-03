%% ADMM method
% A, b: Ax <=b
% opts: optimization parameters
% f: original objective function, sum of areas
% g: original gradient
% h: original hessian 
% r: boundary 

% ======Output========
% x: optimal solution
% opt: objective value
% n: norm of the constraint, stopping criteria


function [x, opt, n] = ADMM(x0, A, b, opts, f, g, h, r, solver, s_opts)
fprintf("ADMM");
fprintf("ITER ; OBJ.VAL; GRAD.NORM");
x = x0;
z = opts.z0;
lambda = opts.lambda0;
iter = 0;
obj_fun = @ (x, z, lamda) f(x) + lambda' * (A*x - z) + opts.sigma/2*norm(A*x - z)^2;
while norm(A*x - z) > opts.tol
    fprintf("iteration start\n")
    
    object = @(x) Object(x, f, z, lambda, opts.sigma, A);
    gradient = @(x) Gradient(x, g, z, lambda, opts.sigma, A);
    hessian = @(x) Hessian(x, h,opts.sigma, A);
    
    % globalized Newton's method
%     newton_opts = struct('maxit', 1e3, 'tol', 1e-5, 'beta1', 1e-6, 'beta2', 1e-6, 'p', 0.1, 'gamma', 0.4, 'sigma', 0.5, 's', 1);
    [x, ~,~, ~] = newton_globalized(object, gradient, hessian, x, s_opts);
    
    %L_BFGS
     % s_opts = struct('s', 1, 'sigma', 0.5, 'gamma', 0.1, 'm', 10, 'epsilon', 1e-7, 'delta', 1e-4,'eta', 0.5, 'H', eye(length(x)));
     % x = solver(object, gradient, x0, s_opts);
%     [ x, ~, ~ ] = PW_L_BFGS (object, gradient, x0, s_opts);
    
    
    z = prox(x, A, b, lambda, opts.sigma);
    lambda = lambda + opts.gamma * opts.sigma * (A* x - z);
    opt = obj_fun(x, z, lambda);
    n = norm(A*addbd(x, r) - z); % vector norm of the constraint
    fprintf("[%4i] ; %2.6f ;%2.6f", iter, n, opt);
    iter = iter + 1;
end
end

function [res] = Object(x, f, z, lambda, sigma, A)
    res = f(x) + lambda' * (A*x - z) + sigma/2*norm(A*x - z)^2;
end

function [res] = Gradient(x, g, z, lambda,sigma, A)
    res = g(x) + sigma* A'*A*x + A'*lambda - A'*z;
end

function [res] = Hessian(x, h,sigma, A)
    res = h(x) + sigma*A'*A;
end


