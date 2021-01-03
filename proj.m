function y = proj(x, A, b)
if A * x <= b
    y = x;
else
    y = x - A' * inv(A*A') * (A * x - b);
end
end