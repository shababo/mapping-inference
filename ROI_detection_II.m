%% ROI detection: Stage II
% Stimulate the selected sites from Stage I

%% set RNG seed
rng(11111,'twister');


% merged_regional
% merged_local
% merged_grid

%% Obtain the location of the new trials 
% 
new_trials = [merged_regional merged_grid(1,3)*ones(size(merged_regional,1),1); merged_grid];
%new_trials = [new_trials; new_trials+normrnd(0,2,size(new_trials)); new_trials+normrnd(0,2,size(new_trials))];
    % merged_grid];
size(new_trials)
%% Prepare and generate new data
% 
num_repeats = 20;

new_trial_locations_on_grid = repmat(new_trials,num_repeats,1);
new_N= size(new_trial_locations_on_grid,1);
pi_k = zeros(new_N,size(all_locations,1));

B = diag(all_locations*inv(A)*all_locations');
for n = 1:new_N
   this_trial_locations = new_trial_locations_on_grid(n,:);
   diff_mat=B+(ones(size(all_locations,1),1)*diag(this_trial_locations*inv(A)*this_trial_locations')')-2*all_locations*inv(A)*this_trial_locations';
    pi_kr = exp(-0.5*(diff_mat));
   pi_k(n,:) = min(.95,sum(pi_kr,2));
end

pi_k_spike = pi_k;
pi_k_spike(pi_k_spike > .65) = 1; % what does this mean?

% % firing delay means and variances
% d_mean_nk = d_mean0 + (1 - pi_nk)*d_mean_coef;
% d_sigma_nk = d_sigma0 + (1 - pi_nk)*d_sigma_coef;
% 
% % sample "ground truth" firing delay
% D = normrnd(d_mean_nk,d_sigma_nk)/data_params.dt + evoked_params.stim_start;
% D(D < evoked_params.stim_start + .002) = evoked_params.stim_start + .002;

% firing delay means and variances
d_mean_nk = d_mean0 + (1.5 - pi_k)*d_mean_coef;
d_sigma_nk = d_sigma0 + (1 - pi_k)*d_sigma_coef;

% sample "ground truth" firing delay
D = exprnd(d_mean_nk/data_params.dt) + evoked_params.stim_start;
% D(D < evoked_params.stim_start + .002) = evoked_params.stim_start + .002;

% sample "ground truth" stimulations

new_X = rand(new_N,K) < pi_k_spike; %.2 
new_X(D > 2000) = 0;

%% Generate a response, given D, Pi, X, w

new_Y = zeros(new_N,data_params.T);

for n = 1:new_N
    firing_neurons = new_X(n,:) & all_amplitudes' > 0;
    
    if any(all_amplitudes(firing_neurons) > 0)
        
        evoked_params.times = D(n,firing_neurons);
        evoked_params.a = all_amplitudes(firing_neurons);
        evoked_params.tau_r = all_tau_rise(firing_neurons)/data_params.dt;
        evoked_params.tau_f = all_tau_fall(firing_neurons)/data_params.dt;
    else
        evoked_params.times = [];
        evoked_params.a = [];
        evoked_params.tau_r = [];
        evoked_params.tau_f = [];
    end
    
    [new_Y(n,:), new_mpp_n] = gen_trace(data_params,bg_params,evoked_params);
    if n == 1
        new_mpp = new_mpp_n;
    else
        new_mpp(n) = new_mpp_n;
    end
end
%% Extract the stimulus-induced events 
new_related_mpp = new_mpp;
new_unrelated_mpp = new_mpp;

for i = 1:new_N
    if size(new_mpp(i).event_times,2) > 0
        indices = new_mpp(i).event_times>evoked_params.stim_start  & new_mpp(i).event_times< (400+evoked_params.stim_start);
        new_related_mpp(i).amplitudes = new_mpp(i).amplitudes(indices);
        new_related_mpp(i).event_times = new_mpp(i).event_times(indices);
        new_unrelated_mpp(i).amplitudes = new_mpp(i).amplitudes(~indices);
        new_unrelated_mpp(i).event_times = new_mpp(i).event_times(~indices);
    end 
end


%% Run analysis on the new data
%% Obtain the background rate of firing, and the distribution of amplitudes 
% 
background_amplitudes = [new_unrelated_mpp.amplitudes];
background_rate = size(background_amplitudes,2)/ (new_N*(data_params.T-400));  % at 100 Hz
background_mean = mean(background_amplitudes);
histogram(background_amplitudes,30)

%% 
n_site = size(new_trials,1);

num_threshold=10;
amplitude_threshold = quantile([new_related_mpp.amplitudes], (1/num_threshold)*[0:num_threshold]);


% Obtain the background rates and averages
background_mean_categories = zeros(num_threshold-1,1);
background_rate_categories = zeros(num_threshold-1,1);
for j = 1:(num_threshold-1)
    amp_selected = background_amplitudes(background_amplitudes>amplitude_threshold(j) & background_amplitudes<(amplitude_threshold(j+2)+0.01));
    background_rate_categories(j) = size(amp_selected,2)/ (new_N*(data_params.T-400));
    background_mean_categories(j) = mean(amp_selected);
end


induced_mean_categories = zeros(n_site,num_threshold-1);
induced_prob_categories = zeros(n_site,num_threshold-1);
induced_prob_categories_sd = zeros(n_site,num_threshold-1);

for n = 1:n_site
    trials_at_this_site = new_related_mpp( n+ ((1:num_repeats)-1)*n_site);
   for j = 1:(num_threshold-1)
         temp = zeros(num_repeats,1);
         temp_sum= zeros(num_repeats,1);
        for i = 1:num_repeats
            selected_amp = trials_at_this_site(i).amplitudes>amplitude_threshold(j)  & trials_at_this_site(i).amplitudes<(amplitude_threshold(j+2)+0.01);
            temp(i) = sum(selected_amp) ;
            temp_sum(i) = sum(  trials_at_this_site(i).amplitudes(selected_amp) ) ;
        end
        bin_index = temp>0;
        prob_index = bin_index -1 + exp(-400*background_rate_categories(j))/exp(-400*background_rate_categories(j));
    induced_prob_categories(n,j) = mean(prob_index);    
    induced_prob_categories_sd(n,j) = std(prob_index)/sqrt(num_repeats);
        if induced_prob_categories_sd(n,j) == 0
            induced_prob_categories_sd(n,j) = sqrt(induced_prob_categories(n,j)*(1-induced_prob_categories(n,j))/num_repeats);
        end
        % To estimate the mean amplitudes, we need to estimate the mean
        % rate first
        mean_overall = mean(temp);
        mean_induced = mean_overall - 400*background_rate_categories(j);
        
        induced_mean_categories(n,j) = ( mean(temp_sum)- background_mean_categories(j)*background_rate_categories(j)*400)/mean_induced;
    end
end



%% When the locations of neurons are known
% 
large_constant = 100;
neuron_likelihood_sum = cell(num_layers,1);
neuron_weighted_amp = cell(num_layers,1);
neuron_likelihood_local_min = cell(num_layers,1);
for i = 1:num_layers
    neuron_likelihood_sum{i}= zeros(size(neuron_locations{i},1),num_threshold -1);
    neuron_likelihood_max{i}=zeros(size(neuron_locations{i},1),num_threshold -1);
end


for j = 1:(num_threshold -1)
   estimates = induced_prob_categories(:,j);
    threshold = 2*induced_prob_categories_sd(:,j);
    % we would need a better standard deviation estimates
    probability = estimates(estimates>threshold);
    prob_sd = induced_prob_categories_sd(estimates>threshold,j);
    locations = new_trials(estimates > threshold,:);
    mean_amp= induced_mean_categories(estimates > threshold,j);
    reference_prob = min(max(probability),1);
    
    % Note here A is the scaling matrix!
    for i = 1:num_layers
        for l = 1:size(neuron_locations{i},1)
            this_one = neuron_locations{i}(l,1:2);
            temp_weights=0;
            temp_weighted_amp = 0;
            for k = 1:size(probability,1)
                dist_scaled = (locations(k,1)-this_one(1))^2/A(1,1)+(locations(k,2)- this_one(2))^2/A(2,2);
                p_ijk = reference_prob*exp(-0.5*dist_scaled);
                % then check the Gaussian density
                
                if prob_sd(k)==0
                    neuron_likelihood_sum{i}(l,j) = neuron_likelihood_sum{i}(l,j)+ large_constant;
                    temp_weights = temp_weights+large_constant;
                    temp_weighted_amp = temp_weighted_amp+mean_amp(k)*large_constant;
                else 
                neuron_likelihood_sum{i}(l,j) = neuron_likelihood_sum{i}(l,j)+normpdf(p_ijk,probability(k),prob_sd(k));
                                   temp_weights = temp_weights+normpdf(p_ijk,probability(k),prob_sd(k));
                    temp_weighted_amp = temp_weighted_amp+mean_amp(k)*normpdf(p_ijk,probability(k),prob_sd(k));
                end
            end
            neuron_weighted_amp{i}(l,j)=temp_weighted_amp/temp_weights;
        end
    end
end

%% Draw the inferred neurons and their synaptic strengths
% Merge the neurons in different layers into one big matrix:

all_likelihood_sum =zeros(0,num_threshold-1);
all_weighted_amp =zeros(0,num_threshold-1);

all_locations =zeros(0,3);
for i = 1:num_layers
    all_likelihood_sum = [all_likelihood_sum; neuron_likelihood_sum{i}];
    all_weighted_amp = [all_weighted_amp; neuron_weighted_amp{i}];
    all_locations = [all_locations; neuron_locations{i}(:,1:3)];
end

all_dist = diag(all_locations(:,1:2)*all_locations(:,1:2)')*ones(1,size(all_locations,1))...
    +ones(size(all_locations,1),1)*diag(all_locations(:,1:2)*all_locations(:,1:2)')'...
    -2*all_locations(:,1:2)*all_locations(:,1:2)';

%% Finding local maximum 
% We let the local maximum be the neurons with largests likelihood among
% its neighbours, given an acceptable radius R
% The distance matrix is calculated in the previous block 
R=2000;
figure(11)
threshold_level=3;
colormap = jet(2);
colormap(1,:)= [1 1 1];
colormap(2,:)= [1 0 0];
for i = 1:num_layers
    connected_neurons_ind = find(neuron_features(i).amplitude);
    temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
        -neuron_locations{i}(connected_neurons_ind,2),...
        neuron_features(i).amplitude(connected_neurons_ind)*25);
    set(temp,'MarkerFaceColor','k');
    alpha(temp,0.6);
    hold on
end
set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
view(2)

    xlim([20,460]);
    ylim([-900,-400]);
    
for j = 1:(num_threshold-1)
neuron_centers=zeros(0,2);
neuron_amps = zeros(0,1);
    likelihood_thresholds =quantile( all_likelihood_sum(:,j), [0.8  0.9 0.95]);
    for i = 1:size(all_dist,1)
        if all_likelihood_sum(i,j) == max(all_likelihood_sum(all_dist(i,:)<R,j))
            if all_likelihood_sum(i,j)>likelihood_thresholds(threshold_level)
            neuron_centers = [neuron_centers; all_locations(i,1:2)];
            neuron_amps = [neuron_amps all_weighted_amp(i,j)];
            
            end
        end
    end
     %potential_neuron = scatter(neuron_centers(:,1),-neuron_centers(:,2),amplitude_threshold(j+1)*15,colormap(2,:),'filled','o');
     
     potential_neuron = scatter(neuron_centers(:,1),-neuron_centers(:,2),neuron_amps*25,colormap(2,:),'filled','o');
     alpha(potential_neuron,0.5);
     
end
hold off
saveas(11,'../Data/Final.jpg')


