%%
%% General L-BFGS Solver
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%% USAGE:
%%  1.       s, initial value of line searching
%%  2.   sigma, line searching exponential scaling factor
%%  3.   gamma, line searching condition scaling factor
%%  4.      x0, initial point
%%  5.       m, memory factor
%%  6. epsilon, terminating tolerance
%%  8.       f, objective function
%%  9.       g, gradient function
%% 10.       H, initial matrix of memoization
%%


function [ opt, x, G ] = L_BFGS (s, sigma, gamma, x0, m, epsilon, f, g, H)
    S    = {}; % store all s^k = difference of x
    Y    = {}; % store all y^k = defierence of gradient
    RHO  = {}; % store all 1/<s^k, y^k> to reduce the computation
    G    = {}; % record gradient norm
    x    = x0; % initialize x
    
    %% main procedure
    gradient   = g(x);
    n          = norm(gradient);
    while n > epsilon
        q          = gradient; 
        G{end + 1} = n;
        [~, l]     = size(S);
        a          = zeros(l, 1);
        for i=l:-1:1
            a(i) = RHO{i} * transpose(S{i}) * q;
            q    = q - a(i) * Y{i};
        end
        r = H * q;
        for i=1:l
            beta = RHO{i} * transpose(Y{i}) * r;
            r = r + (a(i) - beta) * S{i};
        end
        d          = -r; % descent direction
        alpha      = backtracking(s, sigma, gamma, x, d, gradient, f);
        next_x     = x + alpha * d;
        next_s     = next_x - x;
        next_g     = g(next_x);
        next_y     = next_g - gradient;
        
        %% update for next round
        x          = next_x;
        gradient   = next_g; 
        n          = norm(next_g);
        S{end+1}   = next_s;
        Y{end+1}   = next_y;
        RHO{end+1} = 1.0/(transpose(next_s) * next_y);
        [~, l] = size(S);
        if l > m
            S = S(:, 2:end);
            Y = Y(:, 2:end);
            RHO = RHO(:, 2:end);
        end
    end
    opt = f(x);
end