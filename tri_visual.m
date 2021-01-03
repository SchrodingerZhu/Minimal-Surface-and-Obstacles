%%
%% Visualization of Trigularized Surface
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%% USAGE:
%% 1. x0, x1, starting and ending position of x-axis
%% 2. y0, y1, starting and ending position of y-axis
%% 3.      Z, the vector of the display area (dim = m * n)
%%

function [] = tri_visual (x0, x1, y0, y1, Z) 
    [m, n]  = size(Z);
    x       = x0:(x1 - x0)/(m-1):x1;
    y       = y0:(y1 - y0)/(n-1):y1;
    [X, Y]  = meshgrid(x, y);
    T       = delaunay(X, Y);
    trimesh(T, X, Y, Z);
    xlabel("$x$", "interpreter", "latex")
    ylabel("$y$", "interpreter", "latex")
%     colorbar;
end