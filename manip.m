function y = manip(X)
    % Inintialize parameters
    [m, n] = size(X);  % the matrix after adding the boundary. 
    L1 = 1;
    L2 = 1;
    a = L1/(m-1);
    b = L2/(n-1);
    % c = L1*L2 / (m-1)*(n-1);
    c = a*b;
    
    M1 = X(2:end, 1:end-1); % The original matrix
    M2 = X(1:end-1, 1:end-1); % The matrix shifting down
    M3 = X(1:end-1, 2:end); % The matrix shifting left
    M4 = X(2:end, 2:end); % shift up

    s = (m-1) * (n-1); % number of elements in this matrix

    % possible problem: the sign of delta 
    D1 = reshape(M1 - M2, [1, s]); % Delta_down indexed by {2,...,m} and {1,...,n-1}
    D2 = reshape(M1 - M4, [1, s]); % Delta_right indexed by {2,...,m} and {2,...,n} (?)
    D3 = reshape(M3 - M2, [1, s]); % Delta_right indexed by {1,...,m-1} and {2,...,n}
    D4 = reshape(M3 - M4, [1, s]); % Delta_down indexed by {2,...,m} and {2,...,n}
    C = ones(1,s);
    F = [a*D1; b*D2; c*C];
    G = [b*D3; a*D4; c*C];
    y = [F,G];
    % each column corresponds to a vector, we can take norm directly. 
end