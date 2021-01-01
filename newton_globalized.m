% globalized Newton's method
function [x, obj, f_k, grad_k] = newton_globalized(f, grad, hess, x0, opts)

% === INPUT ==========
% f      objective function
% grad   gradient of objective function
% hess   hessian of objective function
% x0     initial point
% opts   a struct for the options:
%          - .maxit         maximum number iterations
%          - .tol           tolerance
%          - .beta1         direction choice parameter
%          - .beta2         direction choice parameter
%          - .p             direction choice parameter
%          - .gamma         line search parameter
%          - .sigma         line search parameter
%          - .s             line search parameter

% === OUTPUT =========
% x      an optimal solution of min f(x)
% obj    the optimal function value up to a tolerance
% f_k    a vector that stores objective function value at each iteration
% grad_k a vector that stores norm of gradient at each iteration
    
    fprintf("− − − globalized Newton's method;\n");
    fprintf("ITER ; OBJ.VAL ; G.NORM ; DIR ; STEP.SIZE\n");
    
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
        
        % stopping criterion
        if ng <= opts.tol
            fprintf("[%4i] ; %1.6f ; %1.6f ; STOP ; %1.4f\n", iter, f_val, ng, alpha);
            break
        end
        
        % compute Newton direction
        s = -hess(x) \ g;
        
        % determine whether or not to use Newton direction
        ns = norm(s);
        if -g' * s >= min(opts.beta1, opts.beta2*ns^opts.p) * ns^2
            d = s;
            fprintf("[%4i] ; %1.6f ; %1.6f ; newton ; %1.4f\n", iter, f_val, ng, alpha);
        else
            d = -g;
            fprintf("[%4i] ; %1.6f ; %1.6f ; gradient ; %1.4f\n", iter, f_val, ng, alpha);
        end
        
        % Armijo backtracking line search
        alpha = backtracking(opts.s, opts.sigma, opts.gamma, x, d, g, f);
        x = x + alpha*d;
    end
    
	obj = f(x);
end