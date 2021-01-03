m = 50;
n = 50;



dim_x = (m-2) * (n-2);

% boundary function
R = { @(x,y) (x^2 - y^2), @(x, y) 1 + sin(2 * pi * x), @(x, y) 1 + cos(1 / (x + 0.001) ), @(x, y) 0.5 - abs(y - 0.5), @(x, y) (1 + exp(xy))^-1, @(x, y) 1 + asin(-1 + 2 * sqrt(x * y)) }


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
saveas(gcf, "obstacles.png")
% initial point

x0 = ones(dim_x, 1);

% constraints
G = @(x)(b - A * x);
GG = @(x)(-transpose(A));
HG = @(x)(sparse(zeros(dim_x)));
H = @(x) 0;
GH = @(x) 0;
HH = @(x) 0;

% PW_L_BFGS
pwl_opts = struct('s', 5, 'sigma', 0.5, 'gamma', 0.1, 'm', 50, 'epsilon', 5e-4, 'eta', 0.5, 'H', eye(dim_x)); 

% globalized Newton's method
newton_opts = struct('maxit', 1e3, 'tol', 1e-5, 'beta1', 1e-6, 'beta2', 1e-6, 'p', 0.1, 'gamma', 0.4, 'sigma', 0.5, 's', 1);

%QPM parameters
penalty = 5;


for i = 1:numel(R)
  r = R{i};

  eps = 1e-6;
  % objective functions
  obj = @(x) objective(x, m, n, r);
  grad = @(x) gradient(x, m, n, r);
  hess = @(x) hessian(x, m, n, r);

  tic
  [x, opt, ng, fs, iter] = qpm_smn(obj, grad, hess, G, GG, H, GH, x0, newton_opts, penalty, eps);
  t2 = toc

  fid = fopen(sprintf("smn%d.txt", i), "wb");
  fwrite(fid, opt);
  fwrite(fid, ng);
  fwrite(fid, fs);
  fwrite(fid, iter);
  fclose(fid);

  figure(i + 6);
  ANS = reshape(x, [m-2, n-2]);
  tri_visual (0, 1, 0, 1, addbd(ANS, r)); 
  saveas(gcf, sprintf("surface-smn%d.png", i));

  figure(i + 7);
  plot(i + ng);
  saveas(gcf, sprintf("norm-smn%d.png", i));

  figure(i + 8);
  plot(fs);
  saveas(gcf, sprintf("obj-smn%d.png", i));

  figure(i + 9);
  plot(iter);
  saveas(gcf, sprintf("iteration-smn%d.png", i));
  

  eps = 1e-3;
  tic
  [x, opt, ng, fs, iter] = qpm(obj, grad, G, GG, H, GH, x0, @PW_L_BFGS, pwl_opts, penalty, eps)
  t1 = toc
  
  fid = fopen(sprintf("bfgs%d.txt", i), "wb");
  fwrite(fid, opt);
  fwrite(fid, ng);
  fwrite(fid, fs);
  fwrite(fid, iter);
  fclose(fid);

  figure(i + 2);
  ANS = reshape(x, [m-2, n-2]);
  tri_visual (0, 1, 0, 1, addbd(ANS, r)); 
  saveas(gcf, sprintf("surface-bfgs%d.png", i));

  figure(i + 3);
  plot(i + ng);
  saveas(gcf, sprintf("norm-bfgs%d.png", i));

  figure(i + 4);
  plot(fs);
  saveas(gcf, sprintf("obj-bfgs%d.png", i));

  figure(i + 5);
  plot(iter);
  saveas(gcf, sprintf("iteration-bfgs%d.png", i));

end

