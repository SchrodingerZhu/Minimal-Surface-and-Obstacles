function y = g(x, m, n, r)
    if length(x) ~= (m-2)*(n-2)
        fprintf("Wrong dimension!\n")
    end
    X = reshape(x, m-2, []);
    y = manip2(addbd(X, r))';
    y = y / 2;
end