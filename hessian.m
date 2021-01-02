% hessian
function H = hessian(x, m, n, r)
    if length(x) ~= (m-2)*(n-2)
        fprintf("Wrong dimension!\n")
    end
    X = reshape(x, m-2, []);
    H = hess_mat(addbd(X, r));
end

% hessian with matrix input
function H = hess_mat(X)
% === INPUT ==========
% X: m-by-n variable with boundary condition

% === OUTPUT ==========
% hess: resulting hessian
    
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
    A0 = sqrt(a^2*(Delta_down(2:m,1:n-1).^2) + b^2*(Delta_right(2:m,2:n).^2) + c^2*ones(m-1,n-1));
    B0 = sqrt(a^2*(Delta_down(2:m,2:n).^2) + b^2*(Delta_right(1:m-1,2:n).^2) + c^2*ones(m-1,n-1));
    
    % augment A and B with zeros to form m-by-n matrices
    % (i,j) entries of A and B are A_(i,j) and B_(i,j) respectively
    A = [zeros(1,n); A0, zeros(m-1,1)];
    B = [zeros(1,n); B0, zeros(m-1,1)];
    
    % matrices used repeatedly (size: (m-2)-by-(n-2))
    RA1 = 1./A(2:m-1, 2:n-1);   % (k,ell) entry is 1/A(k,ell)
    RA2 = 1./A(3:m, 2:n-1);   % (k,ell) entry is 1/A(k+1,ell)
    RA3 = 1./A(2:m-1,1:n-2);   % (k,ell) entry is 1/A(k,ell-1)
    
    RB1 = 1./B(3:m, 1:n-2);   % (k,ell) entry is 1/B(k+1,ell-1)
    RB2 = 1./B(2:m-1, 1:n-2);   % (k,ell) entry is 1/B(k,ell-1)
    RB3 = 1./B(3:m, 2:n-1);   % (k,ell) entry is 1/B(k+1,ell)
    
    X1 = X(2:m-1, 2:n-1);   % (k,ell) entry is X(k,ell)
    X2 = X(1:m-2, 2:n-1);   % (k,ell) entry is X(k-1,ell)
    X3 = X(2:m-1, 3:n);   % (k,ell) entry is X(k,ell+1)
    X4 = X(3:m, 2:n-1);   % (k,ell) entry is X(k+1,ell)
    X5 = X(2:m-1, 1:n-2);   % (k,ell) entry is X(k,ell-1)
    
    % full reciprocal matrix m-by-n
    RA = [zeros(1,n); 1./A0, zeros(m-1,1)];   % (i,j) entry is 1/A(i,j) or 0
    RB = [zeros(1,n); 1./B0, zeros(m-1,1)];   % (i,j) entry is 1/B(i,j) or 0
    
    % six components of the first-order derivative (size: (m-2)-by-(n-2))
    % (i,j) entries are derivatives w.r.t. X(i,j)
    DA1 = 1 ./ A(2:m-1,2:n-1) .* (a^2*Delta_down(2:m-1,2:n-1) - b^2*Delta_right(2:m-1,3:n));
%     DA1 = [zeros(1,n); zeros(m-2,1), DA1, zeros(m-2,1); zeros(1,n)];
    
    DA2 = -a^2 ./ A(3:m,2:n-1) .* Delta_down(3:m,2:n-1);
%     DA2 = [zeros(m-2,1), DA2, zeros(m-2,1); zeros(2,n)];
    
    DA3 = b^2 ./ A(2:m-1,1:n-2) .* Delta_right(2:m-1,2:n-1);
%     DA3 = [zeros(1,n); zeros(m-2,2), DA3; zeros(1,n)];
    
    DB1 = 1 ./ B(3:m,1:n-2) .* (b^2*Delta_right(2:m-1,2:n-1) - a^2*Delta_down(3:m,2:n-1));
%     DB1 = [zeros(m-2,2), DB1; zeros(2,n)];
    
    DB2 = a^2 ./ B(2:m-1,1:n-2) .* Delta_down(2:m-1,2:n-1);
%     DB2 = [zeros(1,n); zeros(m-2,2), DB2; zeros(1,n)];
    
    DB3 = -b^2 ./ B(3:m,2:n-1) .* Delta_right(2:m-1,3:n);
%     DB3 = [zeros(m-2,1), DB3, zeros(m-2,1); zeros(2,n)];
    
    % the negative of derivatives of 1./A and 1./B
    % (i,j) entries are derivatives w.r.t. X(i,j)
    CA1 = (RA1.^2).*DA1;
    CA2 = (RA2.^2).*DA2;
    CA3 = (RA3.^2).*DA3;
    
    CB1 = (RB1.^2).*DB1;
    CB2 = (RB2.^2).*DB2;
    CB3 = (RB3.^2).*DB3;
    
    % step 1: 2nd-order derivative w.r.t. X(k,ell) and X(k,ell)
    S11 = (a^2+b^2)*RA1 + a^2*RA2 + b^2*RA3 + (a^2+b^2)*RB1 + a^2*RB2 + b^2*RB3;
    S12 = -a^2 * X1 .* (CA1 + CA2 + CB1 + CB2) - b^2 * X1 .* (CA1 + CA3 + CB1 + CB3);
    S13 = a^2 * X2 .* (CA1 + CB2);
    S14 = b^2 * X3 .* (CA1 + CB3);
    S15 = a^2 * X4 .* (CA2 + CB1);
    S16 = b^2 * X5 .* (CA3 + CB1);
    
    S1 = (1/2) * (S11 + S12 + S13 + S14 + S15 + S16);
    diag1 = reshape(S1, [], 1);
    
    % add zeros on the boundary for convenience of indexing (size: m-by-n)
    CA1 = [zeros(1,n); zeros(m-2,1), CA1, zeros(m-2,1); zeros(1,n)];
    CA2 = [zeros(1,n); zeros(m-2,1), CA2, zeros(m-2,1); zeros(1,n)];
    CA3 = [zeros(1,n); zeros(m-2,1), CA3, zeros(m-2,1); zeros(1,n)];
    CB1 = [zeros(1,n); zeros(m-2,1), CB1, zeros(m-2,1); zeros(1,n)];
    CB2 = [zeros(1,n); zeros(m-2,1), CB2, zeros(m-2,1); zeros(1,n)];
    CB3 = [zeros(1,n); zeros(m-2,1), CB3, zeros(m-2,1); zeros(1,n)];
    
    % step 2: 2nd-order derivative w.r.t. X(k,ell) and X(k-1,ell)
    S21 = -a^2 * (RA(3:m-1,2:n-1) + RB(3:m-1,1:n-2));
    S22 = -a^2 * Delta_down(3:m-1,2:n-1) .* (CA2(2:m-2,2:n-1) + CB1(2:m-2,2:n-1));
    S23 = b^2 * Delta_right(3:m-1,3:n) .* CA2(2:m-2,2:n-1);
    
    S2 = (1/2) * (S21 + S22 + S23);
    [~, num_col_S2] = size(S2);
    S2 = [S2; zeros(1, num_col_S2)];
    diag2 = reshape(S2, [], 1);
    diag2 = diag2(1:end-1);
%     disp(diag2)
    
    % step 3: 2nd-order derivative w.r.t. X(k,ell) and X(k,ell-1)
    S31 = -b^2 * (RA(2:m-1,2:n-2) + RB(3:m,2:n-2));
    S32 = a^2 * Delta_down(3:m,3:n-1) .* CB3(2:m-1,2:n-2);
    S33 = -b^2 * Delta_right(2:m-1,3:n-1) .* (CA1(2:m-1,2:n-2) + CB3(2:m-1,2:n-2));
    
    S3 = (1/2) * (S31 + S32 + S33);
    diag3 = reshape(S3, [], 1);
%     disp(diag3)
    
    % step 4: 2nd-order derivative w.r.t. X(k,ell) and X(k-1,ell-1)
    S41 = -a^2 * Delta_down(3:m-1,3:n-1) .* CB3(2:m-2,2:n-2);
    S42 = -b^2 * Delta_right(3:m-1,3:n-1) .* CA2(2:m-2,2:n-2);
    
    S4 = (1/2) * (S41 + S42);
%     disp(size(S4))
    [~, num_col_S4] = size(S4);
    S4 = [S4; zeros(1, num_col_S4)];
    diag4 = reshape(S4, [], 1);
    diag4 = diag4(1:end-1);
    
    dim_x = (m-2) * (n-2);
    diag2 = [zeros(dim_x-length(diag2) ,1); diag2];
    diag3 = [zeros(dim_x-length(diag3) ,1); diag3];
    diag4 = [zeros(dim_x-length(diag4) ,1); diag4];
    
    H1 = spdiags(diag1, 0, dim_x, dim_x);
    H2 = spdiags(diag2, 1, dim_x, dim_x);
    H3 = spdiags(diag3, m-2, dim_x, dim_x);
    H4 = spdiags(diag4, m-1, dim_x, dim_x);
%     disp(full(H1))
%     disp(full(H2))
%     disp(full(H3))
    
    H = H1 + H2 + H2' + H3 + H3' + H4 + H4';
end