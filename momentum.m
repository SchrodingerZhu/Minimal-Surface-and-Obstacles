% opts.maxit = 20000; opts.tol = 10^(-5);
% opts.gamma = 0.1; opts.s = 0.001;
% opts.sigma = 0.5; opts.m = 1; opts.e = 10^(-7);
function [x,obj,iter] = momentum(f,x0,opts)
iter = 0;
x1 = [0, x0];
beta = 2/(iter+2); % Other rules to select step size
while norm(f.grad(x1(2)) > opts.tol
    y = x1(2) + beta*(x1(2) - x1(1));
    d = -f.grad(y);
    alpha = opts.s;
    i = 0;
    while (f.obj(y + alpha * d) - f.obj(y) >= -opts.gamma*opts.alpha* (norm(d)^2)% &&(i <10^4)
        alpha = alpha * opts.sigma;
        i = i+1;
    end
%     i
%     alpha
    xtemp = y + alpha*d;
    
    x1(1) = x1(2);
    x1(2) = xtemp;
    
    iter = iter + 1;
end
    
    x = x1(2);
    obj = f.obj(x);
end