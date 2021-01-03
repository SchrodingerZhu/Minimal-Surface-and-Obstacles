% Barzilai-Borwein gradient method with non-monotone line search
function [x, obj, f_k, grad_k, T, method] = gm_bb_nonmonotone(f, grad, x0, opts)

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
% x       an optimal solution of min f(x)
% obj     the optimal function value up to a tolerance
% f_k     a vector that stores objective function value at each iteration
% grad_k  a vector that stores norm of gradient at each iteration
% T       a vector that stores elapsed cpu-time at each iteration 
% method  name of the optimization method
    
    fprintf("− − − Barzalai-Borwein gradient method with nonmonotone line search;\n");
    fprintf("ITER ; OBJ.VAL ; G.NORM ; STEP.SIZE\n");
    method = "Barzalai-Borwein Gradient Method with Nonmonotone Line Search";
    
    f_k = [];
    grad_k = zeros(opts.maxit, 1);
    T = zeros(opts.maxit, 1);  

    % Armijo backtracking line search in the first step
    s = 1;
    sigma = 0.5;
    gamma = 0.5;
    
    tic
    x_old = x0;
    f_old = f(x_old);
    g_old = grad(x_old);
    ng_old = norm(g_old);
    lambda = backtracking(s, sigma, gamma, x_old, -g_old, g_old, f);
    x = x_old - lambda*g_old;
    T(2) = toc;
    
    fprintf("[%4i] ; %2.6f ; %2.6f ; %1.4f\n", 0, f_old, ng_old, lambda);
    
    f_k(end+1) = f_old;
    grad_k(1) = ng_old;
    
    % main loop
    for iter = 2 : opts.maxit
        tic
        g_old = grad(x_old);
        g = grad(x);
        ng = norm(g);
        fprintf("[%4i] ; %2.6f ; %2.6f ; %1.4f\n", iter-1, f(x), ng, lambda);
        
        % store iterates
        f_k(end+1) = f(x);
        grad_k(iter) = ng;
        
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
            alpha = opts.delta;
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
        
        % record time
        T(iter+1) = toc;
    end
    
    obj = f(x);
    f_k = f_k(1 : iter);
    grad_k = grad_k(1 : iter);
    T = T(1 : iter);
    T = cumsum(T);
end