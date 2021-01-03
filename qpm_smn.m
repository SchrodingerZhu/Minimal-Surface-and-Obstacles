%%
%% Quadratic Penalty Method
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

function [x, opt, ng, fs, iter] = qpm_smn(f, g, h, G, GG, HG, H, GH, HH, x0, opts, penalty, eps)
    alpha = penalty;
    p = @(x)(P(x, f, G, H, alpha));
    grad = @(x)(Grad(x, g, G, GG, H, GH, alpha));
    hess = @(x)(Hess(x, h, G, GG, HG, H, GH, HH, alpha, eps))
    x = newton_globalized(p, grad, hess, x0, opts);
    fs = [f(x)];
    ng = [norm(max(0, G(x)))];
    iter = [];
    while not(check(x, G, H, eps))
        alpha = alpha * penalty;
        p = @(x)(P(x, f, G, H, alpha));
        grad = @(x)(Grad(x, g, G, GG, H, GH, alpha));
        hess = @(x)(Hess(x, h, G, GG, HG, H, GH, HH, alpha, eps))
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

function [res] = Hess(x, h, G, GG, HG, H, GH, HH, penalty, eps)
    P1 = h(x);
    % gh = GH(x);
    gg = sparse(GG(x));
    g  = G(x);
    % P2 = gh(:, 1:end) * gh(:, 1:end)'; 
    % P3 = permute(H(x)' .* permute(HH(x), [1 3 2]), [1 3 2]);
    % P4 = permute(max(0, G(x))' .* permute(HG(x), [1 3 2]), [1 3 2]);
    D  = (g <= eps & g >= -eps);
    U1  = gg * diag(g > eps);
    U2  = repmat(D', length(x), 1);
    U   = sparse(U1 + U2);
    P5 = gg(:, 1:end) * U(:, 1:end)';
    res = P1 + penalty * P5; 
end

