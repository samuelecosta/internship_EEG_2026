function G = compute_gcv_matrix(alpha, s2, Y_hat, rho2, m)
    f = alpha ./ (s2 + alpha); 

    num = sum(sum((f .* Y_hat).^2)) + rho2;
    den = (m - sum(s2 ./ (s2 + alpha)))^2;
    G = num / den;
end