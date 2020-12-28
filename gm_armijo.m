% gradient method with Armijo backtracking line search
function [x, obj, f_k, grad_k] = gm_armijo(f, grad, x0,opts)

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
% x      an optimal solution of min f(x)
% obj    the optimal function value up to a tolerance
% f_k    a vector that stores objective function value at each iteration
% grad_k a vector that stores norm of gradient at each iteration
    
    fprintf("− − − gradient method with backtracking;\n");
    fprintf("ITER ; OBJ.VAL ; G.NORM ; STEP.SIZE\n");
    
    % initialization
    x = x0;
    f_k = [];
    grad_k = [];
    alpha = 0;

    % main loop
    for iter = 1 : opts.maxit
        f_val = f(x);
        g = grad(x);
        ng = norm(g);
        fprintf("[%4i] ; %2.6f ; %2.6f ; %1.4f\n", iter-1, f_val, ng, alpha);
        
        % store iterates
        f_k(end+1) = f_val;
        grad_k(end+1) = ng;
        
        % stopping criterion
        if ng <= opts.tol
            break;
        end
        
        % Armijo backtracking line search
        alpha = backtracking(opts.s, opts.sigma, opts.gamma, x, -g, g, f);
        x = x - alpha*g;
    end
    
    obj = f_val;
end