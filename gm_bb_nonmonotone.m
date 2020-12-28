% Barzilai-Borwein gradient method with non-monotone line search
function [x, obj, f_k, grad_k] = gm_bb_nonmonotone(f, grad, x0, opts)

% === INPUT ==========
% f      objective function
% grad   gradient of objective function
% x0     initial point
% opts   a struct for the options:
%          - .maxit     maximum number iterations
%          - .tol       tolerance
%          - .alpha_lb  lower bound for appropriate step size
%          - .alpha_ub  upper bound for appropriate step size
%          - .delta     replacement of inappropriate step size
%          - .gamma     line search parameter
%          - .sigma     line search parameter
%          - .M         number of past iterates used

% === OUTPUT =========
% x      an optimal solution of min f(x)
% obj    the optimal function value up to a tolerance
% f_k    a vector that stores objective function value at each iteration
% grad_k a vector that stores norm of gradient at each iteration
    
    fprintf("− − − Barzalai-Borwein gradient method with nonmonotone line search;\n");
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
    lambda = backtracking(s, sigma, gamma, x_old, -g_old, g_old, f);
    x = x_old - lambda*g_old;
    
    % main loop
    for iter = 2 : opts.maxit
        g_old = grad(x_old);
        g = grad(x);
        ng = norm(g);
        fprintf("[%4i] ; %2.6f ; %2.6f ; %1.4f\n", iter-1, f(x), ng, lambda);
        
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
        alpha = (s'*y) / (s'*s);
        
        % check whether alpha is of appropriate size
        if alpha >= opts.alpha_ub || alpha <= opts.alpha_lb
            alpha = delta;
        end
        
        % nonmonotone line search
        lambda = 1 / alpha;
        f_past = max(f_k(max(1,end-opts.M+1) : end));
        
        factor = -opts.gamma * (g'*g);
        lhs = f(x - lambda*g) - f_past;
        rhs = lambda * factor;
        
        while lhs > rhs
            lambda = lambda * opts.sigma;
            lhs = f(x - lambda*g) - f_past;
            rhs = opts.sigma * rhs;
        end
        
        % update x_old and x
        x_old = x;
        x = x - lambda*g;
    end
    
    obj = f(x);
end