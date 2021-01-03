m = 15;
n = 15;
dim_x = (m-2) * (n-2);

% boundary function
r = @(x,y) x+y;

%% cone obstacle
[X, Y, Z] = get_obstacle(0.8, 0.8, 0.1, 'cone', struct('radius', 0.5, 'height', 6));
scatter3(X, Y, Z);
[A,    b] = obstacles(X, Y, Z, [], [], m - 2, n - 2)

% initial point
x0 = ones(dim_x, 1);

% objective functions
obj = @(x) f(x, m, n, r);
grad = @(x) g(x, m, n, r);


% constraints
% G = @(x)(b - A * x)
% GG = @(x)(-transpose(A))
% H = @(x) 0
% GH = @(x) 0

% PW_L_BFGS
s_opts = struct('s', 1, 'sigma', 0.5, 'gamma', 0.1, 'm', 25, 'epsilon', 1e-4, 'eta', 0.5, 'H', eye(dim_x)); 

% Projected gradient
p_opts = struct('s', 1, 'sigma', 0.5, 'gamma', 0.8, 'tol', 1, 'lambda0', 5, 'L0', 0.1);


[x, opt] = proj_grad(obj, grad, x0, p_opts, A, b);


%QPM parameters
% penalty = 5;
% eps = 1e-3;

% [x, opt] = qpm(obj, grad, G, GG, H, GH, x0, @PW_L_BFGS, s_opts, penalty, eps);

figure(2)
ANS = reshape(x, [m-2, n-2]);
tri_visual (0, 1, 0, 1, addbd(ANS, r)); 
