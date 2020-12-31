function [x,objv,iter] = momentum(obj, grad, x0,opts)
fprintf("momentum");
fprintf("ITER ; OBJ.VAL; STEP.SIZE\n");
iter = 0;
x_old = zeros(size(x0));
x = x0;
% Other rules to select step size
while norm(grad(x)) > opts.tol
    beta = 2/(iter+2);
    %y = x + beta*(x - x_old);
    d = -grad(x);
    alpha = opts.s;
    i = 0;
    while (obj(x + alpha * d) - obj(x) >= -opts.gamma * alpha * norm(d)^2) && (i < 5)
        alpha = alpha * opts.sigma;
        i = i+1;
    end
    xtemp = x + alpha * d + beta*(x - x_old);
    x_old = x;
    x = xtemp;
    iter = iter + 1;
    fprintf("[%4i] ; %2.6f ; %1.4f\n", iter, obj(x),  alpha)
end
    objv = obj(x);
end