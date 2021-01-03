% gradient method with Armijo backtracking line search
function [x, obj, f_k, grad_k, T, method] = gm_armijo(f, grad, x0,opts)

% === INPUT ==========
% f      objective function
% grad   gradient of objective function
% x0     initial point
% opts   a struct for the options:
%          - .maxit         maximum number iterations
%          - .tol           tolerance
%          - .gamma         line search parameter
%          - .sigma         line search parameter
%          - .s             line search parameter

% === OUTPUT =========
% x       an optimal solution of min f(x)
% obj     the optimal function value up to a tolerance
% f_k     a vector that stores objective function value at each iteration
% grad_k  a vector that stores norm of gradient at each iteration
% T       a vector that stores elapsed cpu-time at each iteration 
% method  name of the optimization method
    
    fprintf("− − − gradient method with backtracking;\n");
    fprintf("ITER ; OBJ.VAL ; G.NORM ; STEP.SIZE\n");
    method = "Gradient Method with Armijo Backtracking Line Search";
    
    % initialization
    x = x0;
    alpha = 0;
    f_k = zeros(opts.maxit, 1);
    grad_k = zeros(opts.maxit, 1);
    T = zeros(opts.maxit, 1);

    % main loop
    for iter = 1 : opts.maxit
        tic
        f_val = f(x);
        g = grad(x);
        ng = norm(g);
        fprintf("[%4i] ; %2.6f ; %2.6f ; %1.4f\n", iter-1, f_val, ng, alpha);
        
        % store iterates
        f_k(iter) = f_val;
        grad_k(iter) = ng;
        
        % stopping criterion
        if ng <= opts.tol
            break;
        end
        
        % Armijo backtracking line search
        alpha = backtracking(opts.s, opts.sigma, opts.gamma, x, -g, g, f);
        x = x - alpha*g;
        
        % record time
        T(iter+1) = toc;
    end
    
    obj = f_val;
    f_k = f_k(1 : iter);
    grad_k = grad_k(1 : iter);
    T = T(1 : iter);
    T = cumsum(T);
end