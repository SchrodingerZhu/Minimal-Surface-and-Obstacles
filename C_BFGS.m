%%
%% Compact BFGS Solver
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%% USAGE:
%%  1.       f, objective function
%%  2.       g, gradient function
%%  3.      x0, initial point
%%  4.    opts, options
%%
%% OPTIONS:
%%  1.       s, initial value of line searching
%%  2.   sigma, line searching exponential scaling factor
%%  3.   gamma, line searching condition scaling factor
%%  4.       m, memory factor
%%  5. epsilon, terminating tolerance
%%  6.     eta, pw search lower bound scaling factor
%%

function [ x, opt, G ] = C_BFGS(f, g, x0, opts)
    s       = opts.s;
    sigma   = opts.sigma;
    gamma   = opts.gamma; 
    m       = opts.m;
    epsilon = opts.epsilon;
    eta     = opts.eta;
    % memoization
    
    
    %% initialization with one step of gradient descent
    
    % iteration sensitive quantity
    px       = x0;
    pgrad    = g(x0); 
    d        = -pgrad;
    alpha    = backtracking(s, sigma, gamma, px, d, pgrad, f);
    x        = x0 + alpha * d;
    grad     = g(x);
    n        = norm(grad);
    
    % memoization
    G = {};
    S = {};
    Y = {};
    
    % state transformation
    R        = [];       % previous R
    D        = [];       % previous D
    YY       = []; % previous Y * Y
    
    %% main procedure
    iter = 0;
    opt = f(x);
    fprintf("ITER ; OBJ.VAL ; G.NORM ; STEP.SIZE\n");
    while n > epsilon
        % update memoization
        G{end+1} = n;                   
        S{end+1} = x - px;
        Y{end+1} = grad - pgrad;
        [~, l]   = size(S);
        if l > m 
            S = S(:, 2:end);
            Y = Y(:, 2:end);
        end
        [~, k]   = size(S);
        
        % step one
        gg = transpose(grad) * grad;
        Sg = transpose(cell2mat(S)) * grad;
        Yg = transpose(cell2mat(Y)) * grad;
        
        % update D
        D(end+1, end+1) = transpose(grad - pgrad)*S{end};
        if l > m
            D = D(2:end, 2:end);
        end
        
        % update R
        expR = zeros(k, 1);
        for i=1:k
            expR(i) = transpose(S{i}) * (grad - pgrad);
        end
        if l > m
            R(1:end-1, 1:end-1) = R(2:end, 2:end);
            R(:, end) = expR;
        else
            R(1:k, end+1) = expR;
        end
        
        % update YY
        expYY = zeros(k, 1);
        for i=1:k-1
            expYY(i) = transpose(Y{i}) * (grad - pgrad);
        end
        expYY(k) = -transpose(grad) * grad + 2*transpose(grad - pgrad)*grad + transpose(pgrad)*pgrad;
        if l > m
            YY(1:end-1, 1:end-1) = YY(2:end, 2:end);
            YY(:, end) = expYY;
            YY(end, :) = transpose(expYY);
        else
            YY(1:k, end+1) = expYY;
            YY(end, :) = transpose(expYY);
        end
        
        % update gamma_k
        gamma_k = (transpose(Y{end}) * S{end})/(transpose(Y{end}) * Y{end});
        
        % calculate p and d
        iR = inv(R);
        tR = transpose(R);
        tiR = transpose(iR);
        p = [tiR*(D + gamma_k*YY)*iR*Sg - gamma_k * tiR * Yg; -iR*Sg];
        d = -gamma_k*grad - [cell2mat(S), gamma_k * cell2mat(Y)] * p;
        
        % backtracking and update
        alpha = pw_search(s, sigma, gamma, eta, x, d, f, g, epsilon);
        pgrad = grad;
        px    = x;
        x     = x + alpha * d;
        grad  = g(x);
        n     = norm(grad);
        iter  = iter + 1;
        opt = f(x);
        fprintf("[%4i] ; %2.6f ; %2.6f ; %1.4f\n", iter, opt, n, alpha);
    end
    
    
end