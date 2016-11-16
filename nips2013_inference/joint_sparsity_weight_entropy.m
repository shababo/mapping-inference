function entropy = joint_sparsity_weight_entropy(alpha, mu, sig_sq)

log_alpha = log(alpha);
log_alpha_less_1 = log(1-alpha);

log_alpha(isinf(log_alpha)) = 0;
log_alpha_less_1(isinf(log_alpha_less_1)) = 0;

entropy_sparsity = -sum(alpha.*log_alpha + (1 - alpha).*log_alpha_less_1);

entropy_weights = sum(alpha .* (0.5 * (1 + log(2*pi*sig_sq))));

entropy = entropy_weights + entropy_sparsity;
