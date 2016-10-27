function H_expected = approximate_expected_joint_entropy_multiple_neurons(x_next,a, m, rho_squared, sigma_n, options, num_samples)

if nargin < 6
    num_samples = 100;
end

K = length(a);

H_expected = 0;

for sample = 1:num_samples
    % sample y_next
    gamma_sample = rand(K,1) < a;
    w_sample = gamma_sample.*normrnd(m,sqrt(rho_squared));
    y = normrnd(x_next*w_sample,sigma_n);
    
    % update params
    [alpha_sample, mu_sample, s_sq_sample] = run_varbvs_general(x_next, y, sigma_n, rho_squared, a, m, options);
    
    % calculate per-sample entropy
    H_expected = H_expected + sum(per_neuron_joint_entropy(alpha_sample, mu_sample, s_sq_sample));
end

% compute expectation
H_expected = H_expected/num_samples;