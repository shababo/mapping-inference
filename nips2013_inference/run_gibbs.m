function output = run_gibbs(params,data)

a = .05*ones(params.K,1);

c = ones(params.K,1);

K = params.K;
N = params.N;
%% Generate some data

% these parameters govern the time delay, as a function of the
% point-spread-function stimuli for a particular trial
d_mean0 = 65;
d_sigma0 = 5;
d_mean_coef = 50;
d_sigma_coef = 20;


% probability of firing
pi_kr = exp(-0.5*squareform(pdist(params.coords,'mahalanobis',params.A)).^2);
pi_nk = zeros(params.N,params.K);
for n = 1:params.N
    pi_nk(n,:) = min(1,sum(pi_kr.*data.stims(n*ones(params.K, 1),:),2)');
end

% firing delay means and variances
d_mean_nk = d_mean0 + (1 - pi_nk)*d_mean_coef;
d_sigma_nk = d_sigma0 + (1 - pi_nk)*d_sigma_coef;



%% now we can gibbs sample it and see how we do. initialize sampler:

% L samples
L = 100;

% B burn-in samples
B = 40;

% store post-burnin samples here:
w_samples = zeros(K,L);
gamma_samples = zeros(K,L);
X_samples = zeros(N,K,L);
D_samples = zeros(N,K,L);

% current iteration / initialization
gamma_s = rand(K,1) < a;
w_s = gamma_s .* (sign(c - .5) .* abs(normrnd(0,1,K,1)));
D_s = d_mean0 + (1 - pi_nk)*d_mean_coef;
X_s = rand(N,K) < pi_nk;


%% okay, let's sample

for sample = -B:L

    fprintf('Sample %d of %d\n', sample, L);
    [X_s, D_s, w_s, gamma_s] = gibbs_single_sweep(X_s, D_s, w_s, gamma_s, data.responses, pi_nk, c, a, params.sigma_s, params.sigma_n, d_mean_nk, d_sigma_nk, params.t, params.tau, params.g);

    if sample > 0
        
        w_samples(:,sample) = w_s;
        gamma_samples(:,sample) = gamma_s;
        D_samples(:,:,sample) = D_s;
        X_samples(:,:,sample) = X_s;
        
        figure(4); bar(mean(w_samples(:,1:sample),2));
        drawnow;
        
    end

end

output.w_samples = w_samples;
output.gamma_samples = gamma_samples;
output.D_samples = D_samples;
output.X_samples = X_samples;



