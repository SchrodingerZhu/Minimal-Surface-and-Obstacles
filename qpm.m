%%
%% Quadratic Penalty Method
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%% USAGE:
%%  1.       f, objective function
%%  2.       g, gradient of objective function 
%%  3.       G, inequality constraint functions  (cell array)
%%  4.      GG, inequality gradients (cell array)    
%%  5.       H, equality constraint functions (cell array)
%%  6.      GH, equality gradients (cell array)
%%  7.      x0, initial point
%%  8.  solver, unconstraint solver
%%  9.  s_opts, solver options
%% 10. penalty, penalty parameter (> 1)
%% 11,     eps, constraint checking tolerance
%%

function [x, opt] = qpm(f, g, G, GG, H, GH, x0, solver, s_opts, penalty, eps)
    alpha = penalty;
    p = @(x)(P(x, f, g, G, H, alpha));
    grad = @(x)(Grad(x, g, G, GG, H, GH, alpha));
    x = solver(p, grad, x0, s_opts);
    while not(check(x, G, H, eps))
        alpha = alpha * penalty;
        p = @(x)(P(x, f, g, G, H, alpha));
        grad = @(x)(Grad(x, g, G, GG, H, GH, alpha));
        x = solver(p, grad, x0, s_opts);
    end
    opt = f(x);
end

function [res] = check(x, G, H, eps)
    res = true;
    for i = 1:length(G)
        res = res && (G{i}(x) <= eps);
    end
    for i = 1:length(H)
        res = res && (H{i}(x) <= eps) && (H{i}(x) >= -eps);
    end
end

function [res] = P(x, f, G, H, penalty)
    res = f(x);
    for i = 1:length(G)
        res = res + 0.5 * penalty * max(0, G{i}(x))^2;
    end
    for i = 1:length(H)
        res = res + 0.5 * penalty * H{i}(x)^2;
    end
end

function [res] = Grad(x, g, G, GG, H, GH, penalty)
    res = g(x);
    for i = 1:length(GG)
        res = res + penalty * max(0, G{i}(x)) * GG{i}(x);
    end
    for i = 1:length(H)
        res = res + penalty * H{i}(x) * GH{i}(x);
    end
end

