% augmented Lagrange method for problems with equality constraints
function [x, obj] = augment_lagrange(f, h, grad_f, grad_h, x0, mu0, solver, opts)
% === INPUT ==========
% f       objective function
% h       equality constraints h(x)=0
% grad_f  gradient of f
% grad_h  gradients of the equalities, column i as grad h_i
% x0      initial point
% mu0     initial multiplier
% solver  solver in each inner loop
% opts    a struct for the options:
%          - .maxit     maximum number iterations
%          - .tol       tolerance, a function of the iteration number
%          - .alpha     penalty parameter, a function of the iteration number

% === OUTPUT =========
% x      an optimal solution of min f(x)
% obj    the optimal function value up to a tolerance
% f_k    a vector that stores objective function value at each iteration

    fprintf("− − − augmented Lagrange method;\n");
    fprintf("ITER ; OBJ.VAL\n");
    
    % initialization
    x = x0;
    mu = mu0;
    obj_k = f(x);
    
    % main loop
    for k = 1 : opts.maxit
        fprintf("[%4i] ; %2.6f\n", k, obj_k)
        
        % options for inner minimization of the augmented Lagrange function
        inner_opts.tol = opts.tol(k);
        inner_opts.maxit = 1e3;
        inner_opts.gamma = 0.5;
        inner_opts.sigma = 0.5;
        inner_opts.s = 1;
        
        % minimize augmented Lagrange function
        aug_L = @(x) f(x) + mu'*h(x) + opts.alpha(k)/2 * (h(x)'*h(x));
        grad_aug_L = @(x) grad_f(x) + grad_h(x)*mu + opts.alpha(k)*grad_h(x)*h(x);
        [x, obj_k, ~, ~] = gm_armijo(aug_L, grad_aug_L, x, inner_opts);
        
        % ascent step for multiplier
        mu = mu + opts.alpha(k)*h(x);
    end
    
    obj = f(x);
end