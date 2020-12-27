%%
%% Backtracking Linear Searcher
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%% USAGE:
%% 1.     s, initial value of line searching
%% 2. sigma, line searching exponential scaling factor
%% 3. gamma, line searching condition scaling factor
%% 4.     x, current position
%% 5.     d, descent direction
%% 6.  grad, gradient vector at x
%% 7.     f, target function
%%
%% OUTPUT:
%% alpha, a scalar represents the step size
%%

function [alpha] = backtracking(s, sigma, gamma, x, d, grad, f) 
    alpha = s;
    factor = gamma * transpose(d) * grad;
    value = f(x);
    lhs = f(x + alpha * d) - value;
    rhs = alpha * factor;
    while lhs > rhs
        alpha = alpha * sigma;
        lhs = f(x + alpha * d) - value;
        rhs = alpha * factor;
    end
end