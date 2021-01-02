% objective function
function y = objective(x, m, n, r)
    if length(x) ~= (m-2)*(n-2)
        fprintf("Wrong dimension!\n")
    end
    X = reshape(x, m-2, []);
    y = sum(f_mat(addbd(X, r)));
end

% objective function with matrix input
function obj = f_mat(X)
% === INPUT ==========
% X: m-by-n variable with boundary condition

% === OUTPUT ==========
% f: resulting objective function
    
    [m, n] = size(X);
    L1 = 1;
    L2 = 1;
    a = L1 / (n-1);
    b = L2 / (m-1);
    c = a * b;
    s = (m-1) * (n-1);
    
    X_right = [zeros(m,1), X(:,1:end-1)];
    X_down = [zeros(1, n); X(1:m-1,:)];

    Delta_right = X - X_right;
    Delta_down = X - X_down;
    C = c * ones(s,1);
    
    M1 = [a * reshape(Delta_down(2:m,1:n-1), [s,1]), b * reshape(Delta_right(2:m,2:n), [s,1]), C];
    M2 = [a * reshape(Delta_down(2:m,2:n), [s,1]), b * reshape(Delta_right(1:m-1,2:n), [s,1]), C];
    
    M1 = M1';
    M2 = M2';
    
    obj = sum(0.5 * sqrt(sum([M1,M2].^2)));   % vector norm on each column
end