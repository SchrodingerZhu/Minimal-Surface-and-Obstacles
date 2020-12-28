% pure Barzilai-Borwein gradient method
function [x, obj, f_k, grad_k] = gm_bb(f, grad, x0, opts)

% === INPUT ==========
% f      objective function
% grad   gradient of objective function
% x0     initial point
% opts   a struct for the options:
%          - .maxit  maximum number iterations
%          - .tol    tolerance
%          - .step   choice of step size (1 or 2)

% === OUTPUT =========
% x      an optimal solution of min f(x)
% obj    the optimal function value up to a tolerance
% f_k    a vector that stores objective function value at each iteration
% grad_k a vector that stores norm of gradient at each iteration
    
    fprintf("− − − Barzalai-Borwein gradient method with step size choice %u;\n", opts.step);
    fprintf("ITER ; OBJ.VAL ; G.NORM ; STEP.SIZE\n");
    
    f_k = [];
    grad_k = [];
    
    % initialization
    x_old = x0;
    f_old = f(x_old);
    g_old = grad(x_old);
    ng_old = norm(g_old);
    
    fprintf("[%4i] ; %2.6f ; %2.6f\n", 0, f_old, ng_old);
    f_k(end+1) = f_old;
    grad_k(end+1) = ng_old;

    % Armijo backtracking line search in the first step
    s = 1;
    sigma = 0.5;
    gamma = 0.5;
    alpha = backtracking(s, sigma, gamma, x_old, -g_old, g_old, f);
    x = x_old - alpha*g_old;
    
    % main loop
    for iter = 2 : opts.maxit
        g_old = grad(x_old);
        g = grad(x);
        ng = norm(g);
        fprintf("[%4i] ; %2.6f ; %2.6f ; %1.4f\n", iter-1, f(x), ng, alpha);
        
        % store iterates
        f_k(end+1) = f(x);
        grad_k(end+1) = ng;
        
        % stopping criterion
        if ng <= opts.tol
            break;
        end

        % Barzilai-Borwein step size
        s = x - x_old;
        y = g - g_old;
        
        if opts.step == 1
            alpha = (s'*s) / (s'*y);
        else
            alpha = (s'*y) / (y'*y);
        end
        
        % update x_old and x
        x_old = x;
        x = x - alpha*g;
    end
    
    obj = f(x);
end