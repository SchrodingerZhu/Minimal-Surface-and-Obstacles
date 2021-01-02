% gradient
function grad = gradient(x, m, n, r)
    if length(x) ~= (m-2)*(n-2)
        fprintf("Wrong dimension!\n")
    end
    X = reshape(x, m-2, []);
    grad = g_mat(addbd(X, r));
end

% gradient with matrix input
function g = g_mat(X)
% === INPUT ==========
% X: m-by-n variable with boundary condition

% === OUTPUT ==========
% g: resulting gradient
    
    [m, n] = size(X);
    L1 = 1;
    L2 = 1;
    a = L1 / (n-1);
    b = L2 / (m-1);
    c = a * b;
    
    X_right = [zeros(m,1), X(:,1:end-1)];
    X_down = [zeros(1, n); X(1:m-1,:)];

    Delta_right = X - X_right;
    Delta_down = X - X_down;
    
    % sizes of A and B are both (m-1)-by-(n-1)
    A = sqrt(a^2*(Delta_down(2:m,1:n-1).^2) + b^2*(Delta_right(2:m,2:n).^2) + c^2*ones(m-1,n-1));
    B = sqrt(a^2*(Delta_down(2:m,2:n).^2) + b^2*(Delta_right(1:m-1,2:n).^2) + c^2*ones(m-1,n-1));
    
    % augment A and B with zeros to form m-by-n matrices
    % (i,j) entries of A and B are A_(i,j) and B_(i,j) respectively
    A = [zeros(1,n); A, zeros(m-1,1)];
    B = [zeros(1,n); B, zeros(m-1,1)];
    
    % six components of the gradient
    S1 = 1./A(2:m-1,2:n-1) .* (a^2*Delta_down(2:m-1,2:n-1) - b^2*Delta_right(2:m-1,3:n));
    S2 = b^2 ./ A(2:m-1,1:n-2) .* Delta_right(2:m-1,2:n-1);
    S3 = -a^2 ./ A(3:m,2:n-1) .* Delta_down(3:m,2:n-1);
    
    S4 = 1./B(3:m,1:n-2) .* (b^2*Delta_right(2:m-1,2:n-1) - a^2*Delta_down(3:m,2:n-1));
    S5 = -b^2 ./ B(3:m,2:n-1) .* Delta_right(2:m-1,3:n);
    S6 = a^2 ./ B(2:m-1,1:n-2) .* Delta_down(2:m-1,2:n-1);
    
    S = 1/2 * (S1 + S2 + S3 + S4 + S5 + S6);
    g = reshape(S, [(m-2)*(n-2), 1]);
end