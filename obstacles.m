%%
%% Add Obstable to Current (A, b)
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%%  USAGE:
%%  1.    X, Y, obstacle coordinate
%%  2.       Z, obstacle hight
%%  3.       A, current chooser
%%  4.       b, current limit
%%  5.    m, n, dimensions
%%

function [A, b] = obstacles(X, Y, Z, A, b, m, n)
    X = transpose(X);
    Y = transpose(Y);
    [height, width] = size(A);
    col = n * (X-1) + Y;
    row = height+1:height+length(X);
    A   = transpose([A; sparse(length(X), m * n)]);
    indices = (m * n)*(row - 1) + col;
    A(indices) = 1;
    b = [b; Z];
    A   = transpose(A);
end