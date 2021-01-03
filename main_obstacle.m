m = 50;
n = 50;
dim_x = (m-2) * (n-2);

% boundary function
r = @(x,y) -2 * sin(2 * pi * x) + 2 * sin(2 * pi * y);


VX = [];
VY = [];
VZ = [];
A  = [];
b  = [];
% semishpere obstacle
[X, Y, Z] = get_obstacle(1, 1, 0.1, 'semisphere', struct('radius', 0.7));
[A,    b] = obstacles(X, Y, Z, A, b, m - 2, n - 2);
VX = [VX; X];
VY = [VY; Y];
VZ = [VZ; Z];
% cone obstacle
[X, Y, Z] = get_obstacle(2.8, 2.8, 0.1, 'cone', struct('radius', 0.5, 'height', 6));
[A,    b] = obstacles(X, Y, Z, A, b, m - 2, n - 2);
VX = [VX; X];
VY = [VY; Y];
VZ = [VZ; Z];
%rect obstacle
[X, Y, Z] = get_obstacle(1, 2.8, 0.1, 'rect', struct('length', 0.5, 'width', 0.5, 'height', 6));
[A,    b] = obstacles(X, Y, Z, A, b, m - 2, n - 2);
VX = [VX; X];
VY = [VY; Y];
VZ = [VZ; Z];

figure(1);
scatter3(VX, VY, VZ);
% initial point

x0 = rand(dim_x, 1);

% objective functions
obj = @(x) objective(x, m, n, r);
grad = @(x) gradient(x, m, n, r);
hess = @(x) hessian(x, m, n, r);

% constraints
G = @(x)(b - A * x);
GG = @(x)(-transpose(A));
HG = @(x)(sparse(zeros(dim_x)));
H = @(x) 0;
GH = @(x) 0;
HH = @(x) 0;

% PW_L_BFGS
s_opts = struct('s', 5, 'sigma', 0.5, 'gamma', 0.1, 'm', 25, 'epsilon', 1e-4, 'eta', 0.5, 'H', eye(dim_x)); 

% gradient method with Barzilai-Borwein step and nonmonotone line search
bbn_opts = struct('maxit', 5e3, 'tol', 1e-5, 'alpha_lb', 1e-4, 'alpha_ub', 1e4, 'delta', 1, 'gamma', 0.5, 'sigma', 0.5, 'M', 5);

% globalized Newton's method
newton_opts = struct('maxit', 1e3, 'tol', 1e-5, 'beta1', 1e-6, 'beta2', 1e-6, 'p', 0.1, 'gamma', 0.4, 'sigma', 0.5, 's', 1);

%QPM parameters
penalty = 5;
eps = 1e-3;

%tic
%[x, opt, ng, fs, iter] = qpm(obj, grad, G, GG, H, GH, x0, @PW_L_BFGS, s_opts, penalty, eps)
%toc

tic
[x, opt, ng, fs, iter] = qpm_smn(obj, grad, hess, G, GG, H, GH, x0, newton_opts, penalty, eps)
toc

figure(2);
ANS = reshape(x, [m-2, n-2]);
tri_visual (0, 1, 0, 1, addbd(ANS, r)); 

figure(3);
plot(ng);

figure(4);
plot(fs);

figure(5);
plot(iter);

