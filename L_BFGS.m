%%
%% General L-BFGS Solver
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%% USAGE:
%%  1.       f, objective function
%%  2.       g, gradient function
%%  3.      x0, initial point
%%  4.    opts, options

%% OPTIONS
%%  1.       s, initial value of line searching
%%  2.   sigma, line searching exponential scaling factor
%%  3.   gamma, line searching condition scaling factor
%%  4.       m, memory factor
%%  5. epsilon, terminating tolerance
%%  6.   delta, skipping tolerance
%%  7.       H, initial matrix of memoization
%%


function [ x, opt, F, G, T, method ] = L_BFGS (f, g, x0, opts)
    s       = opts.s;
    sigma   = opts.sigma;
    gamma   = opts.gamma; 
    m       = opts.m;
    epsilon = opts.epsilon;
    delta   = opts.delta;
    H       = opts.H;

    
    S      = {}; % store all s^k = difference of x
    Y      = {}; % store all y^k = difference of gradient
    RHO    = {}; % store all 1/<s^k, y^k> to reduce the computation
    F      = []; % record objective function value
    G      = []; % record gradient norm
    T      = [0]; % record elapsed time
    x      = x0; % initialize x
    method = "L-BFGS Method"; % name of the optimization method
    
    %% main procedure
    gradient   = g(x);
    n          = norm(gradient);
    opt        = f(x);
    iter       = 0;
    fprintf("− − − L-BFGS method;\n");
    fprintf("ITER ; OBJ.VAL ; G.NORM ; STEP.SIZE; SKIP \n");
    while n > epsilon
        tic
        q        = gradient; 
        G(end+1) = n;
        F(end+1) = opt;
        [~, l]   = size(S);
        a        = zeros(l, 1);
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
        
        % update for next round
        x          = next_x;
        gradient   = next_g; 
        n          = norm(next_g);
        rho_inv    = transpose(next_s) * next_y;
        
        % check skipping condition
        if rho_inv >= delta
          S{end+1}   = next_s;
          Y{end+1}   = next_y;
          RHO{end+1} = 1.0/rho_inv;
          [~, l] = size(S);
          if l > m
              S = S(:, 2:end);
              Y = Y(:, 2:end);
              RHO = RHO(:, 2:end);
          end
        end
        
        iter = iter + 1;
        opt = f(x);
        fprintf("[%4i] ; %2.6f ; %2.6f ; %1.4f; %d\n", iter, opt, n, alpha, rho_inv < delta);
        
        T(end+1) = toc;
    end
    
    G(end+1) = n;
    F(end+1) = opt;
    T = cumsum(T);
end
