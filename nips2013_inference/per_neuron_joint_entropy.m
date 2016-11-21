function H_k = per_neuron_joint_entropy(alpha, mu, s_squared)

log_alpha = log(alpha);
log_alpha_less_1 = log(1-alpha);

log_alpha(isinf(log_alpha)) = 0;
log_alpha_less_1(isinf(log_alpha_less_1)) = 0;

entropy_sparsity = -(alpha.*log_alpha + (1 - alpha).*log_alpha_less_1);

entropy_weights = alpha .* (0.5 * (1 + log(2*pi*s_squared)));

H_k = entropy_weights + entropy_sparsity;