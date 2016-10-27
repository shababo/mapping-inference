function H_k = approximate_expected_joint_entropy(a, m, rho_squared, sigma_n, num_samples)

if nargin < 5
    num_samples = 100;
end

K = length(a);

H_k = zeros(num_samples, K);

for sample = 1:num_samples
    % sample y_next
    gamma_sample = rand(K,1) < a;
    y_next = gamma_sample .* normrnd(m, sqrt(sigma_n^2 + rho_squared)) ...
             + (1 - gamma_sample) .* normrnd(0, sigma_n);
    
    % update params
    s_k_squared = (rho_squared.*sigma_n.^2)./(rho_squared + sigma_n.^2);
    mu_k = (m .* sigma_n.^2)./(rho_squared + sigma_n.^2) + (s_k_squared .* y_next)./ sigma_n.^2;
    log_odds_k = a./(1-a) .* sqrt(s_k_squared./rho_squared) .* exp(mu_k.^2./(2*s_k_squared) - m.^2./(2*rho_squared));
    alpha_k = log_odds_k ./ (1 + log_odds_k);

    % deal with +-infinite log odds
    alpha_k(log_odds_k == Inf) = 1;
    alpha_k(log_odds_k == -Inf) = 0;
    
    % calculate per-sample entropy
    H_k(sample, :) = per_neuron_joint_entropy(alpha_k, mu_k, s_k_squared);
end

% compute expectation
H_k = mean(H_k);
H_k = H_k(:);