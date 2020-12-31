function y = f(x, m, n, r) %make the function calculation efficient
    if length(x) ~= (m-2)*(n-2)
        fprintf("Wrong dimension!\n")
    end
    X = reshape(x, m-2, []);
    y = sum(vecnorm(manip(addbd(X, r))));
    y = y/2;
end






