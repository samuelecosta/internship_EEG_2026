function alpha_opt = get_global_alpha(U, s2, Y)
    [m, ~] = size(U);
    Y_hat = U' * Y;

    rho2 = norm(Y, 'fro')^2 - norm(Y_hat, 'fro')^2; 
    if rho2 < 0, rho2 = 0; end
    
    gcv_cost = @(alpha) compute_gcv_matrix(alpha, s2, Y_hat, rho2, m);
    
    alpha_min = min(s2(s2 > 0)) * 1e-6; 
    alpha_max = max(s2); 
    
    options = optimset('TolX', 1e-8, 'Display', 'off');
    alpha_opt = fminbnd(gcv_cost, alpha_min, alpha_max, options);
end