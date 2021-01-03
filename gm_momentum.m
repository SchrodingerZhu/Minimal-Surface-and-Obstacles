% inertial gradient method
function [x, obj, f_k, grad_k, T, method] = gm_momentum(f, grad, x0, opts)
% === INPUT ==========
% f      objective function
% grad   gradient of objective function
% x0     initial point
% opts   a struct for the options:
%          - .beta     momentum parameter
%          - .tol      tolerance
%          - .L0       initial Lipschitz constant
%          - .sigma    line search parameter for Lipschitz constant

% === OUTPUT =========
% x       an optimal solution of min f(x)
% obj     the optimal function value up to a tolerance
% f_k     a vector that stores objective function value at each iteration
% grad_k  a vector that stores norm of gradient at each iteration
% T       a vector that stores elapsed cpu-time at each iteration 
% method  name of the optimization method

    fprintf("− − − inertial gradient method;\n");
    fprintf("ITER ; OBJ.VAL; GRAD.NORM; STEP.SIZE\n");
    method = "Inertial Gradient Method";
    
    iter = 0;
    x_old = zeros(size(x0));
    x = x0;
    L = opts.L0; % initial Liptistz constant
    
    f_k = [];
    grad_k = [];
    T = [0];

    while norm(grad(x)) > opts.tol
        tic
        f_k(end+1) = f(x);
        grad_k(end+1) = norm(grad(x));
        
        y = x + opts.beta*(x - x_old);
        d = grad(x);
        alpha = 2*(1-opts.beta)/L;
        delta = opts.beta*(x - x_old) - alpha * d;
        
        % search for Lipschitz constant
        while f(y - alpha*d) > f(x) + d'*delta + L/2*(norm(delta))^2
            L = L/opts.sigma;
            alpha = 2*(1-opts.beta)/L;
            delta = opts.beta*(x - x_old) - alpha * d;
        end
        
        xtemp = y - alpha * d;
        x_old = x;
        x = xtemp;
        iter = iter + 1;
        fprintf("[%4i] ; %2.6f ;%2.6f; %1.4f\n", iter, f(x), norm(grad(x)), alpha)
        
        T(end+1) = toc;
    end
    
    obj = f(x);
    f_k(end+1) = obj;
    grad_k(end+1) = norm(grad(x));
    T = cumsum(T);
end
