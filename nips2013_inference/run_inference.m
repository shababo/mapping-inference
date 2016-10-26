function output = run_inference(params, data)


pi_kr = exp(-0.5*squareform(pdist(params.coords,'mahalanobis',params.A)).^2);

alpha_sum = sum(alpha_synapse(params.t,0,params.tau,-params.g));

hyperparam_sigma_n = sqrt(length(params.t))*params.sigma_n/abs(alpha_sum);


hyperparam_p_connected = .1*ones(params.K,1);


Y_n = sum(data.responses,2)/alpha_sum;


pi_nk = zeros(params.N,params.K);
for n = 1:params.N
    pi_nk(n,:) = min(1,sum(pi_kr.*data.stims(n*ones(params.K, 1),:),2)');
end

alpha = zeros(params.K,1); %ones(params.K, 1) * alpha_0;
mu = zeros(params.K, 1);
s_sq = zeros(params.K,1);

n_varbvs_samples = 20;

% run_varbvs(X, Y, sigma_n, sigma_s, alpha, options)
% run_varbvs_general(X, Y, sigma_n, sigma_s, alpha, eta, options);
for sample = 1:n_varbvs_samples
    [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs(pi_nk>rand(params.N,params.K), Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1));%, params.eta);
    alpha = alpha+alpha_tmp/n_varbvs_samples;
    mu = mu+mu_tmp/n_varbvs_samples;
    s_sq = s_sq+s_sq_tmp/n_varbvs_samples;
end

output.alpha = alpha;
output.mu = mu;
output.s_sq = mu;

output.pi_kr = pi_kr;
output.pi_nk = pi_nk;
output.Y_scalar = Y_n;


output.w_estimate = alpha.*mu;

