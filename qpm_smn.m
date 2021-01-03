%%
%% Specialized 
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%% USAGE:
%%  1.       f, objective function
%%  2.       g, gradient of objective function 
%%  3.       G, inequality constraint functions  (vector valued)
%%  4.      GG, inequality gradients (vector valued)    
%%  5.       H, equality constraint functions (vector valued)
%%  6.      GH, equality gradients (vector valued)
%%  7.      x0, initial point
%%  8.  solver, unconstraint solver
%%  9.  s_opts, solver options
%% 10. penalty, penalty parameter (> 1)
%% 11,     eps, constraint checking tolerance
%%

function [x, opt, ng, fs, iter] = qpm_smn(f, g, h, G, GG, H, GH, x0, opts, penalty, eps)
    alpha = penalty;
    p = @(x)(P(x, f, G, H, alpha));
    grad = @(x)(Grad(x, g, G, GG, H, GH, alpha));
    hess = @(x)(Hess(x, h, G, GG, alpha, eps));
    x = newton_globalized(p, grad, hess, x0, opts);
    fs = [f(x)];
    ng = [norm(max(0, G(x)))];
    iter = [];
    while not(check(x, G, H, eps))
        alpha = alpha * penalty;
        p = @(x)(P(x, f, G, H, alpha));
        grad = @(x)(Grad(x, g, G, GG, H, GH, alpha));
        hess = @(x)(Hess(x, h, G, GG, alpha, eps));
        [x, ~, r] = newton_globalized(p, grad, hess, x, opts);
        fs(end+1) = f(x);
        ng(end+1) = norm(max(0, G(x)));
        iter(end+1) = numel(r);
    end
    opt = f(x);
end

function [res] = check(x, G, H, eps)
    res = all((G(x) <= eps) & (H(x) <= eps) & (H(x) >= -eps));
end

function [res] = P(x, f, G, H, penalty)
    res = f(x) + 0.5 * penalty * sum(max(0, G(x)).^2) + 0.5 * penalty * sum(H(x).^2);
end

function [res] = Grad(x, g, G, GG, H, GH, penalty)
    res = g(x) + penalty * GG(x) * max(0, G(x)) + penalty * GH(x) * H(x);
end

function [res] = Hess(x, h, G, GG, penalty, eps)
    P1 = h(x);
    gg = sparse(GG(x));
    g  = G(x);
    U1  = sparse(gg * diag(g > 0));
    P5 = gg(:, 1:end) * U1(:, 1:end)' + eps * speye(length(x));
    res = P1 + penalty * P5; 
end

