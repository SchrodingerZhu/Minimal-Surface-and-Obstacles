function [x,objv,iter] = momentum(obj, grad, x0,opts)
fprintf("momentum");
fprintf("ITER ; OBJ.VAL; GRAD.NORM; STEP.SIZE\n");
iter = 0;
x_old = zeros(size(x0));
x = x0;
L = opts.L0; % adjust the liptistz constant
% Other rules to select step size
while norm(grad(x)) > opts.tol
    y = x + opts.beta*(x - x_old);
    d = grad(x);
    alpha = (1-opts.beta)/L;
    i = 0;
    delta = opts.beta*(x - x_old) - alpha * d;
    while obj(y - alpha*d) > obj(x) + d'*delta + L/2*(norm(delta))^2
        L = L/opts.sigma;
        alpha = (1-opts.beta)/L;
        delta = opts.beta*(x - x_old) - alpha * d;
    end
    xtemp = y - alpha * d;
    x_old = x;
    x = xtemp;
    iter = iter + 1;
    fprintf("[%4i] ; %2.6f ;%2.6f; %1.4f\n", iter, obj(x), norm(grad(x)), alpha)
end
    objv = obj(x);
end
