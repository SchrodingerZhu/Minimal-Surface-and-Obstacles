%% Proximal operator, project x onto an indicator function
% phi(x) = 0    , x<=b;
%          infty, o/w
% x -- projected point
% lambda -- multiplier
% sigma -- proximal weight parameter

function y = prox(x, A, b, lambda, sigma)
y = min (A*x + lambda/sigma,  b); % Q: should this min be taken compoment wise?
    
