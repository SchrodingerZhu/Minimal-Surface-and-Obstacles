function D = manip(X)
%global m n L1 L2
% Inintialize parameters
m = 3;
n = 3;
L1 = 1;
L2 = 1;
a = L1/m;
b = L2/n;
c = L1*L2/m*n; %What's the correct formula???
M1 = X(2:end,1:end-1) % The original matrix
M2 = X(1:end-1,1:end-1); % The matrix shifting down
M3 = X(1:end-1,2:end); % The matrix shifting left
s = numel(M1) % number of elements in this matrix

A = reshape(M1 - M2, [1, s]);
B = reshape(M1 - M3, [1, s]);
C = ones(1,s);
D = [b.* ones(1,s).*A ; a.*ones(1, s).*B; c.*ones(1, s).*C]; % each colume
% corresponds to one vector, we can take norm directly. 
