% plot objective function gap and norm of gradient 
% w.r.t number of iterations and elapsed cpu-time
function [] = plot_iter(f_k, grad_k, T, method)

% === INPUT ==========
% f_k     a vector that stores objective function value at each iteration
% grad_k  a vector that stores norm of gradient at each iteration
% T       a vector that stores elapsed cpu-time at each iteration 
% method  name of the optimization method

    num_iter = length(f_k);
    iters = 0 : num_iter - 1;
    obj = f_k(end);
    gap_k = abs(f_k-obj) / max(1, obj);
    
    figure()
%     plot(iters, gap_k, "Color", "#EDB120")
    plot(iters, gap_k, "Marker", ".", "MarkerFaceColor", "#EDB120", "Color", "#EDB120")
    title({"Plot of $|f(X^k)-f^*|/\max\{1,\,|f^*|\}$ v.s. No. of Iterations", method}, "interpreter", "latex")
    xlabel("Iteration Number $k$", "Interpreter", "latex")
    ylabel("$|f(X^k)-f^*|/\max\{1,\,|f^*|\}$", "Interpreter", "latex")
    
    figure()
%     plot(T, gap_k, "Color", "#EDB120")
    plot(T, gap_k, "Marker", ".", "MarkerFaceColor", "#EDB120", "Color", "#EDB120")
    title({"Plot of $|f(X^k)-f^*|/\max\{1,\,|f^*|\}$ v.s. Time", method}, "interpreter", "latex")
    xlabel("Elapsed CPU Time (seconds)", "Interpreter", "latex")
    ylabel("$|f(X^k)-f^*|/\max\{1,\,|f^*|\}$", "Interpreter", "latex")
    
    figure()
%     plot(iters, grad_k, "Color", "#7E2F8E")
    plot(iters, grad_k, "Marker", ".", "MarkerFaceColor", "#7E2F8E", "Color", "#7E2F8E")
    title({"Plot of $\Vert\nabla f(X^k)\Vert$ v.s. No. Iterations", method}, "interpreter", "latex")
    xlabel("Iteration Number $k$", "Interpreter", "latex")
    ylabel("$\Vert\nabla f(X^k)\Vert$", "Interpreter", "latex")
    
    figure()
%     plot(T, grad_k, "Color", "#7E2F8E")
    plot(T, grad_k, "Marker", ".", "MarkerFaceColor", "#7E2F8E", "Color", "#7E2F8E")
    title({"Plot of $\Vert\nabla f(X^k)\Vert$ v.s. Time", method}, "interpreter", "latex")
    xlabel("Elapsed CPU Time (seconds)", "Interpreter", "latex")
    ylabel("$\Vert\nabla f(X^k)\Vert$", "Interpreter", "latex")
