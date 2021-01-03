%% Projected gradient method
% obj -- objective function 
% grad -- gradient
% x0 -- initial point
% opts -- optimization parameters
% A, b -- linear constraint Ax <= b

function [x, opt] = proj_grad(obj, grad, x0, opts, A, b)
fprintf("Projected gradient")
fprintf("ITER ; OBJ.VAL ; PROJ.G.NORM ; STEP.SIZE\n");
iter = 0;
L = opts.L0;

% different step size rules
%lambda = 1/((iter+1)*L);
lambda = 1/(10*L);

x = x0; % a column vector
opt = obj(x);
gv = grad(x);
d = (proj(x - lambda * gv, A, b) - x); % projected gradient direction
% d = - gv;
armd = gv' * d;
nd = norm(d);
 % want step size to be large, move faster, may not converge. 
% Need further investigation

while (nd > lambda * opts.tol) % && (iter <= 100)
     % Armijo line search to determine step size
    alpha = opts.s;
    i = 0;
    while (obj(x + alpha * d) - opt > opts.gamma * alpha * gv' * d) && (i<=12)
        alpha = alpha * opts.sigma;
        i = i + 1;
    end
    xtemp = x + alpha * d;
    
    % Approximate the L constant
    while obj(xtemp) > obj(x) + grad(x)' * alpha * d + L/2 * alpha^2 * norm(d)
        L = L/opts.sigma;
    end
    % step size rules
    %lambda = 1/((iter+1)*L);
    lambda = 1/(10*L);
    x = xtemp;
    opt = obj(x);
    gv = grad(x);
    % lambda = 1/(iter * norm(gv));
    d = (proj(x - lambda*gv, A, b) - x); % projected gradient direction
    % d = - lambda * gv;
    armd = gv' * d
    nd = norm(d);
    fprintf("[%4i]; %2.4f; %2.4f ; %1.4f\n", iter, opt, nd, alpha);
    i;
    iter = iter + 1;
end
end