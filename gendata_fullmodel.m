%% Generate data from circuit mapping model

% This first pass is based on code we used to generate similar data from
% Shababo, Paige et al 2013 and also from code used to generate
% voltage-clamp data from Merel, Shababo et al 2016.

% set RNG seed
rng(12,'twister');

%% Gen neurons and their types/locations/features

%% build layers

% set priors on layer boundaries (from Lefort et al 2009)
num_layers = 7;
layer_names = {'L1','L2','L3','L4','L5A','L5B','L6'};
layer_bottom_means = [0 128 269 418 588 708 890 1154];
layer_bottom_sds = [0 18 36 44 44 62 80 116]/10;

% draw layer boundaries
layer_boundaries = normrnd(layer_bottom_means,layer_bottom_sds);
while any(diff(layer_boundaries) < 20)
    layer_boundaries = normrnd(layer_bottom_means,layer_bottom_sds);
end

% plot layers
% figure(1234)
% plot([0 1],-bsxfun(@times,[ones(num_layers + 1,2)],[layer_boundaries]'))

%% draw number of cells per layer
% set priors (from Lefort et al 2009)
exc_neurons_per_layer_mean = [0 546 1145 1656 454 641 1288];
exc_neurons_per_layer_sd = [0 120 323 203 112 122 205];
K_layers = ceil(normrnd(exc_neurons_per_layer_mean,exc_neurons_per_layer_sd));
while any(K_layers < 0)
    K_layers = ceil(normrnd(exc_neurons_per_layer_mean,exc_neurons_per_layer_sd));
end

% size of region containing neurons (or region we can stim)
barrel_width = 300;
slide_width = 300;

% how many neurons are excitatory (for now we will only consider a
% homogenous population of excitatory neurons - in the future we may
% consider many different cell types whose properties are different
pct_excitatory = 1.00;
neuron_locations = cell(num_layers,1);
for i = 1:num_layers
    neuron_locations{i} = sample_neuron_positions(K_layers(i), ...
        [0 barrel_width; layer_boundaries(i) layer_boundaries(i+1); 0 slide_width]);
end

figure(12341)
for i = 1:num_layers
    scatter3(neuron_locations{i}(:,1),-neuron_locations{i}(:,2),neuron_locations{i}(:,3),'.');
    hold on
end
hold off
%%


% the K-vector "a" holds the base connectivity to the post-synaptic neuron.

% parameters governing probability of neurons being connected, given the
% are excitatory / inhibitory, and the distance between the neurons
% a_excite_lambda = .005;
% a_inhibit_var = 5000;
% prob_a_max_excite = .22;
% prob_a_max_inhibit = .5;
% 
% a = compute_prior_connectivity_probability(c, p, p_star, ...
%        prob_a_max_excite, prob_a_max_inhibit, a_excite_lambda, a_inhibit_var);

a = .015;


% draw latent weights and sparsity pattern

% prior weight variances. we take these as known, at the moment
slab_sd_excite = 20;
slab_sd_inhib = .75;
%slab_mode_excite = 1;
%slab_mode_inhib = -1.5;
slab_mode_excite = 25;
slab_mode_inhib = 0;

% draw the sparsity pattern
gamma = rand(K,1) < a;

% draw weights
w = gamma.*c.*(slab_mode_excite + rnd_truncated_normal(-slab_mode_excite*ones(K,1))*slab_sd_excite);
w(w < 0) = -w(w<0);
% w = w - gamma.*~c.*(-slab_mode_inhib + rnd_truncated_normal(slab_mode_inhib*ones(K,1))*slab_sd_inhib);

% put together slab priors
sigma_s = slab_sd_excite*ones(K,1).*c + slab_sd_inhib*ones(K,1).*~c;

% draw taus for each neuron

%tau priors
%excite
tau_r_bounds = [1 20]/20000;
tau_f_bounds = [50 200]/20000;
%inhib
% tau_r_bounds = [40 150]/20000;
% tau_f_bounds = [200 500]/20000;

% draw taus for each synapse
tau_r_k = unifrnd(tau_r_bounds(1),tau_r_bounds(2),[K,1]);
tau_f_k = unifrnd(tau_f_bounds(1),tau_f_bounds(2),[K,1]);



%% Plot the neurons

figure(1);
scatter3(p(c,2),-p(c,1),p(c,3), 'b.'); hold on;
% scatter3(p(~c,1),p(~c,2),p(~c,3), 'g.'); hold on;
scatter3(0,0,0,'c.'); hold on;
scatter3(p(gamma,2),-p(gamma,1),p(gamma,3),'ro');
hold off;
legend({'excitatory', 'postsynaptic','connected'})
drawnow;
%input('hit enter to continue')

%% voltage-clamp parameters (daq stuff, bg psc parameters)


data_params.T = 2000; %bins - start not too long
data_params.dt = 1/20000; %
data_params.baseline = 0;
data_params.sigmasq = 3.5;
data_params.phi = [1, .80, -.12]; %this determines what the AR noise looks like.


bg_params.tau_r_bounds = tau_r_bounds*20000;
bg_params.tau_f_bounds = tau_f_bounds*20000;
bg_params.a_min = .5;
bg_params.a_max = 10;
bg_params.firing_rate = 20; %spike/sec 





%% stim paramters

% covariance of point spread function
A = diag([225, 225, 10000]);

evoked_params.stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
evoked_params.stim_tau_fall = .013*20000;
evoked_params.stim_amp = 0;
evoked_params.stim_start = .005*20000;
evoked_params.stim_duration = .005*20000;

% effect on postsynaptic cell
evoked_params.stim_tau_rise = .0015*20000; % values for chr2 from lin et al 2009 (biophysics)
evoked_params.stim_tau_fall = .013*20000;
evoked_params.stim_amp = 0;




%% Generate some data

N = 1323;
R = 1;

% these parameters govern the time delay, as a function of the
% point-spread-function stimuli for a particular trial
% in seconds
d_mean0 = .000;
d_sigma0 = .002;
d_mean_coef = .01;
d_sigma_coef = .005;

% the neurons being stimulated at each trial
% Z = false(N,K);
Z = zeros(N/3,3);
trial_grid_locations = zeros(N/3,2);
count = 1;
for i = 1:21
    for j = 1:21
        trial_grid_locations(count,:) = [i j];
        count = count + 1;
        Z((i-1)*21 + j,1) = (i-1)*20 - 200;
        Z((i-1)*21 + j,2) = (j-1)*20 - 200;
        
    end
end
Z = repmat(Z,3,1);
trial_grid_locations = repmat(trial_grid_locations,3,1);

% probability of firing
% pi_nk = zeros(N,K);
pi_kr = exp(-0.5*squareform(pdist([Z; p],'mahalanobis',A)).^2);
pi_nk = pi_kr(1:N,N+1:N+K);
pi_nk(pi_nk > .65) = 1;
% for n = 1:N
%     stimulus = randsample(K,R);
%     Z(n,stimulus) = 1;
%     pi_nk(n,:) = min(1,sum(pi_kr.*Z(n*ones(K, 1),:),2)');
% end

% firing delay means and variances
d_mean_nk = d_mean0 + (1 - pi_nk)*d_mean_coef;
d_sigma_nk = d_sigma0 + (1 - pi_nk)*d_sigma_coef;

% sample "ground truth" firing delay
D = normrnd(d_mean_nk,d_sigma_nk)/data_params.dt + evoked_params.stim_start;
% D = D;

% D(D < 4) = 4;

% sample "ground truth" stimulations

X = rand(N,K) < pi_nk; %.2 
X(D > 2000) = 0;

%% Generate a response, given D, Pi, X, w

% 
% t = 0:1:2000;
% T = length(t);
% Y = zeros(N,T);
% gmax = 1;
% tau = 2.5;
% sigma_n = 2.5;
% 
% for n = 1:N
%     for k = 1:K
%         if gamma(k) == 1 && X(n,k) == 1
%             Y(n,:) = Y(n,:) + w(k)*X(n,k)*alpha_synapse(t,D(n,k),tau,-gmax);
%         end 
%     end
%     Y(n,:) = Y(n,:) + normrnd(0,sigma_n,1,length(t));
% end

% PUT IN NEW CODE FOR EVENTS
    % - NEED TO DRAW TAUS FOR EACH CELL AS WELL

Y = zeros(N,data_params.T);

for n = 1:N
    firing_neurons = X(n,:) & w' > 0;
    
    if any(w(firing_neurons) > 0)
        
        evoked_params.times = D(n,firing_neurons);
        evoked_params.a = w(firing_neurons);
        evoked_params.tau_r = tau_r_k(firing_neurons)/data_params.dt;
        evoked_params.tau_f = tau_f_k(firing_neurons)/data_params.dt;
    else
        evoked_params.times = [];
        evoked_params.a = [];
        evoked_params.tau_r = [];
        evoked_params.tau_f = [];
    end
    
    Y(n,:) = gen_trace(data_params,bg_params,evoked_params);
    
%     if sum(firing_neurons) && sum(evoked_params.a)
%         break
%     end
    
end


%% plot response

Y_grid = unstack_traces(Y,trial_grid_locations);

figure(2)
plot_trace_stack_grid(Y_grid,3,1,0);


%% now we can gibbs sample it and see how we do. initialize sampler:

% L samples
% L = 100;
% 
% % B burn-in samples
% B = -1;
% 
% % store post-burnin samples here:
% w_samples = zeros(K,L);
% gamma_samples = zeros(K,L);
% X_samples = zeros(N,K,L);
% D_samples = zeros(N,K,L);
% 
% % current iteration / initialization
% gamma_s = rand(K,1) < a;
% w_s = gamma_s .* (sign(c - .5) .* abs(normrnd(0,1,K,1)));
% D_s = d_mean0 + (1 - pi_nk)*d_mean_coef;
% X_s = rand(N,K) < pi_nk;
% 
% 
% %% okay, let's sample
% 
% for sample = -B:L
% 
%     fprintf('Sample %d of %d\n', sample, L);
%     [X_s, D_s, w_s, gamma_s] = gibbs_single_sweep(X_s, D_s, w_s, gamma_s, Y, pi_nk, c, a, sigma_s, sigma_n, d_mean_nk, d_sigma_nk, t, tau, gmax);
% 
%     if sample > 0
%         
%         w_samples(:,sample) = w_s;
%         gamma_samples(:,sample) = gamma_s;
%         D_samples(:,:,sample) = D_s;
%         X_samples(:,:,sample) = X_s;
%         
%     end
% %%
%     figure(3); imagesc(X_s); 
%     figure(5); imagesc(D_s); 
%     figure(4); bar([w, mean(w_samples(:,1:sample),2)]); legend({'weights', 'current sample'});
%     title(num2str(norm(w - mean(w_samples(:,1:sample),2))/norm(w)))
%     drawnow;
% end
% 
figure(2)% %%
% figure(6); bar([w, mean(w_samples,2)]); legend({'weights', 'estimate'});


