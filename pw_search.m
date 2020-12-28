%%
%% PW Line Searcher
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%%  USAGE:
%%  1.     s, initial value of line searching
%%  2. sigma, line searching exponential scaling factor
%%  3. gamma, line searching condition scaling factor
%%  4.   eta, lower bound scaling factor
%%  5.     x, current position
%%  6.     d, descent direction
%%  8.     f, target function
%%  9.     g, gradient function
%% 10.   eps, binary search error limit
%%
function [res] = pw_search(s, sigma, gamma, eta, x, d, f, g, eps)
  grad  = g(x);
  l     = backtracking(s, sigma, gamma, x, d, grad, f);
  r     = l / sigma;
  fval  = f(x);
  gd    = transpose(d) * grad;
  check = @(alpha) (pw_condition(gamma, eta, gd, x, d, alpha, fval, f, g));
  res   = binary_search(l, r, check, eps);
end

function [res] = pw_condition(gamma, eta, gd, x, d, alpha, fval, f, g)
    if f(x + alpha * d) - fval > gamma * alpha * gd
        res = 1;
    else
        if transpose(d) * g(x + alpha * d) < eta * gd
            res = -1;
        else 
            res = 0;
        end
    end
    
end


