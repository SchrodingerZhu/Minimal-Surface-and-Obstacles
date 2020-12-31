function G = manip2(X)

    % Inintialize parameters
    [m, n] = size(X);
    L1 = 1;
    L2 = 1;
    b = L1/(m-1);
    a = L2/(n-1);
    %c = L1*L2/(m-1)*(n-1);
    c = a*b;

    M1 = X(2:end,1:end-1); % The original matrix
    M2 = X(1:end-1,1:end-1); % The matrix shifting down
    M3 = X(1:end-1,2:end); % The matrix shifting left
    M4 = X(2: end, 2:end); % shift up

    s = numel(M1); % number of elements in this matrix
    [h,k] = size(M1); % the rows and cols in M1. 

    D1 = M1 - M2;
    D2 = M1 - M4;
    D3 = M3 - M2;
    D4 = M3 - M4;

    D10 = reshape(D1, [1, s]); %delta 1
    D20 = reshape(D2, [1, s]); %delta 2
    D30 = reshape(D3, [1, s]); %delta 3  
    D40 = reshape(D4, [1, s]); %delta 4
    C = ones(1,s);

    A0 = [a*D10 ; b*D20; c*C]; % each colume corresponds to the index difference at one point. 
    B0 = [b*D30; a*D40; c*C];

    A = reshape(vecnorm(A0),[h, k]); %?? What's the dimension before?
    % each entry is the vector norm at the grid point. 
    B = reshape(vecnorm(B0), [h, k]);

    SD1 = B(1:end-1,1:end-1);
    SN1 = -D4(1:end-1, 1:end-1);
    S1 = a^2.*SN1./SD1; % Some problem here? Why is every entry the same?

    SD2 = A(1:end-1, 1:end-1);
    SN2 = -D2(1:end-1, end-1);
    S2 = b^2.*SN2./SD2;

    SD3 = B(2:end, 1:end - 1);
    SN31 = D3(2:end, 1:end - 1);
    SN32 = D4(2:end, 1:end - 1);
    S3 = (b^2.*SN31 + a^2.*SN32)./SD3;

    SD4 = A(2:end, 2:end);
    SN4 = -D1(2:end, 2:end);
    S4 = a^2.*SN4./SD4;

    SD5 = B(2:end, 2:end);
    SN5 = -D3(2:end, 2:end);
    S5 = b^2.*SN5./SD5;

    SD6 = A(1:end-1, 2:end);
    SN61 = D1(1:end-1, 2:end);
    SN62 = D2(1:end-1, 2:end);
    S6 = (a^2.*SN61 + b^2.*SN62)./SD6;

    l = numel(SD1);

    G = reshape(S1 + S2 + S3 + S4 + S5 + S6, [1, l]);
end