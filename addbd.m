% add boundary to the variable
function Y = addbd(X, r)

% === INPUT ==========
% X: variable of dimension (m-2)-by-(n-2)
% r: boundary function r(x,y)

% === OUTPUT ==========
% Y: variable of dimension m-by-n formed by augmenting X with boundary
    
    [m, n] = size(X);
    m = m + 2;
    n = n + 2;
    % corners of the gird
    a1 = 0;
    b1 = 1;
    a2 = 0;
    b2 = 1;
    
    % length and width of the grid
    L1 = b1 - a1;
    L2 = b2 - a2;
    
    % x- and y-coordinates of the nodes
    t1 = a1: L1/(n-1): b1;   % x-coordinate
    t2 = b2-L2/(m-1): -L2/(m-1): a2+L2/(m-1);   % y-coordinate
    
    % boundary
    r_left = @(y) r(a1, y);
    r_right = @(y) r(b1, y);
    r_up = @(x) r(x, b2);
    r_down = @(x) r(x, a2);
    
    bd_left = arrayfun(r_left, t2)';   % left
    bd_right = arrayfun(r_right, t2)';   % right
    bd_up = arrayfun(r_up, t1);   % up
    bd_down = arrayfun(r_down, t1);   % down
    
    % add boundary
    Z = [bd_left, X, bd_right];
    Y = [bd_up; Z; bd_down];
end