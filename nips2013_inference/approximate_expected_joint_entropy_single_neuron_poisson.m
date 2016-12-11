function H_k = approximate_expected_joint_entropy_single_neuron_poisson(a, m, rho_squared, sigma_n, num_samples)

if nargin < 5
    num_samples = 100;
end

m(m<0)=0;

K = length(a);

H_k = zeros(num_samples, K);

for sample = 1:num_samples
    % sample y_next
    gamma_sample = rand(K,1) < a;
    y_next = gamma_sample .* poissrnd(m);
    
    % update params
    s_k_squared = (rho_squared.*sigma_n.^2)./(rho_squared + sigma_n.^2);
    mu_k = (m .* sigma_n.^2)./(rho_squared + sigma_n.^2) + (s_k_squared .* y_next)./ sigma_n.^2;
    odds_k = a./(1-a) .* sqrt(s_k_squared./rho_squared) .* exp(mu_k.^2./(2*s_k_squared) - m.^2./(2*rho_squared));
    alpha_k = odds_k ./ (1 + odds_k);

    % deal with +-infinite log odds
    alpha_k(odds_k == Inf) = 1;
    alpha_k(odds_k == -Inf) = 0;
    
    % calculate per-sample entropy
    H_k(sample, :) = per_neuron_joint_entropy(alpha_k, mu_k, s_k_squared);
end

% compute expectation
H_k = mean(H_k);
H_k = H_k(:);