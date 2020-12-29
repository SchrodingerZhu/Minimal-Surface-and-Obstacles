function Y = addbd(X, m, n)
% Need to rewrite the parameter part
a1 = 0;
b1 = 1;
a2 = 0;
b2 = 1;
L1 = b1 - a1;
L2 = b2 - a2;
t1 = a1:L1/m:b1;
t2 = a2+L2/n:L2/n:b2-L2/n;
bd1 = r(a2.*ones(size(t1)), t1); % left
bd2 = r(b2.*ones(size(t1)), t1); % right
bd3 = r(a1.*ones(size(t2)), t2); % up
bd4 = r(b1.*ones(size(t2)), t2); % down
Z = [bd3; X; bd4];
Y = [bd1, Z, bd2];