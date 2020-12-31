% opts.maxit = 20000; opts.tol = 10^(-5);
% opts.gamma = 0.1; opts.s = 0.001;
% opts.sigma = 0.5; opts.m = 1; opts.e = 10^(-7);
function [x,obj,iter] = momentum(f, grad, x0, opts)
iter = 0;
x = x0;
x_old = zeros(size(x0));
g = grad(x);
ng = norm(g);

fprintf("ITER ; OBJ.VAL ; GRAD.NORM ; STEP.SIZE\n")
while ng > opts.tol
    beta = 2/(iter+2); % Other rules to select step size
    y = x + beta*(x-x_old);
    g = grad(x);
    d = -g;
    ng = norm(g);
    alpha = opts.s;
    i = 0;
    
    while (f(x + alpha * d)) - f(x) >= -opts.gamma*alpha*ng^2% &&(i <10^4)
        alpha = alpha * opts.sigma;
        i = i + 1;
    end
    
    fprintf("NUM.BACKTRACK = %u\n", i)
    
    fprintf("[%4i] ; %2.6f ; %2.6f ; %2.6f\n", iter, f(x), ng, alpha);

    xtemp = y + alpha*d;
    
    x_old = x;
    x = xtemp;
    
    iter = iter + 1;
end
    
    obj = f.obj(x);
end