addpath(genpath('../../psc-detection'),genpath('../../mapping-inference'),genpath('../../mapping-core'));
%% Load the data set
load('./Environments/6_3_s2c3_mpp_and_stim_data.mat')
cell_locations=cell_locs;
n_cell = size(cell_locations,1);
%% Put the cells into ten non-overlapping groups by their z-coordinates
n_planes = 10;
z_quantiles = quantile(cell_locations(:,3), (1:(n_planes))*0.1);
cell_group_idx = zeros(n_cell,1);
for i_cell = 1:size(cell_locations,1)
    cell_group_idx(i_cell)= sum(cell_locations(i_cell,3)>z_quantiles)+1;
end
cell_group_list = cell(n_planes,1);
for i_plane = 1:n_planes
    cell_group_list{i_plane} = find(cell_group_idx==i_plane);
end
%% Consider the second plane in this analysis 
this_plane = 2;
r1=5;r2=10;num_per_grid=4;
grid_jitters = zeros(num_per_grid,2);
for i_grid = 1:num_per_grid
   grid_jitters(i_grid,:)=[sin(2*pi*i_grid/num_per_grid) cos(2*pi*i_grid/num_per_grid)]; 
end
grid_jitters=[grid_jitters zeros(num_per_grid,1)];
% Calculate the stimulation locations 
cell_group_list{this_plane};

target_locations{this_plane} = zeros(length(cell_group_list{this_plane})*(2*num_per_grid+1),3);


for i_cell_index=1:length(cell_group_list{this_plane})
    i_cell= cell_group_list{this_plane}(i_cell_index);
    nucleus_loc=cell_locations(i_cell,:);
    grid_locs=nucleus_loc;
    grid_locs=[grid_locs;nucleus_loc+grid_jitters*r1];
    grid_locs=[grid_locs;nucleus_loc+grid_jitters*r2];
    target_idx=(i_cell_index-1)*(2*num_per_grid+1) +(1: (2*num_per_grid+1));
    target_locations{this_plane}(target_idx,:) = grid_locs;
end
target_locations{this_plane}(:,3)= z_quantiles(this_plane-1)/2+z_quantiles(this_plane)/2;

plot(target_locations{this_plane}(:,2),target_locations{this_plane}(:,1),'.')
%% Calculate the size of stimuli 
% Load the shape template
load('./Environments/l23_template_cell.mat');
load('./Environments/l23_cells_for_sim.mat');
%l23_average_shape
temp=l23_average_shape;temp_max = max(max(max(temp)));
l23_average_shape = temp/temp_max;
shape_template=l23_average_shape;

cell_params.locations =  cell_locations(cell_group_list{this_plane},:);
cell_params.shape_gain = ones(length(cell_group_list{this_plane}),1);
cell_template = struct();
cell_template.shape= shape_template;
% [pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
%     cell_template,target_locations);
[pi_target, inner_normalized_products] = get_weights_v2(cell_params, ...
    cell_template,target_locations{this_plane});

%% Estimate the first spike probability given a template cell:
% Load the current template
load('./Environments/chrome-template-3ms.mat');
downsamp=1;max_time=300;
power_level = [50 75 100];
num_power_level=length(power_level);
current_template=template(1:downsamp:max_time);
t_vect= 1:1:max_time;
% Delay parameters
delay_params.type=2; %1: normal; 2: gamma
delay_params.mean=58; delay_params.std=15;
delay_params.delayed=true; delay_params.n_grid=200;

g=0.02;
stim_unique = (1:1000)/10;
num_stim_unique= length(stim_unique);
v_trace=cell([num_stim_unique 1]);
n_grid=length(current_template);
for i_stim = 1:num_stim_unique
    temp_trace = zeros([n_grid 1]);
    stims_temp=current_template*stim_unique(i_stim);
    temp_trace(1)=0;
    for i_t = 2:n_grid
        temp1=reshape(stims_temp(1:(i_t-1)), [i_t-1,1]);
        temp2=reshape(exp( ((1:(i_t-1)) -i_t+1)*g), [1,i_t-1]);
        temp_trace(i_t) = temp2*temp1;
    end
    v_trace{i_stim}=temp_trace;
    fprintf('%d voltage grid done;\n', i_stim);
end
linkfunc = {@link_sig, @derlink_sig, @invlink_sig,@derinvlink_sig};

delay_prob = zeros(delay_params.n_grid+1,1);
if delay_params.type == 1 % normal
    delay_prob = normpdf((0:delay_params.n_grid),delay_params.mean,...
        delay_params.std);
elseif delay_params.type == 2 % gamma
    shape=(delay_params.mean^2)/(delay_params.std^2);
    %scale
    scale = delay_params.mean/shape;
    delay_prob = gampdf((0:delay_params.n_grid),shape,scale);
end
% we approximate the probability with densities
delay_prob = delay_prob/sum(delay_prob);
min_delay = 0;
max_delay = delay_params.n_grid;
% Evaluate the probability given the gain_grid
prob_trace=zeros([num_stim_unique 1]);
v_th_known_one=15;
gain_template = 0.02;
for i_stim = 1:num_stim_unique
    v_trace_one=v_trace{i_stim};
        gain_one=gain_template;
        [prob_first_spike_delayed] = voltage_to_prob(gain_one,  v_trace_one,...
            v_th_known_one,linkfunc,delay_params,min_delay,max_delay,delay_prob);
        prob_trace(i_stim) =sum(prob_first_spike_delayed);
end

%% Calculate the probability of spiking induced by each location at two power levels 
prob_cells_by_locations=zeros(size(pi_target,2),size(pi_target,1));
for l = 1:size(pi_target,2)
    for i_cell_index = 1:length(cell_group_list{this_plane})
        prob_cells_by_locations(l,i_cell_index)= stim_to_prob(...
            pi_target(i_cell_index,l)*power_level(3),stim_unique,prob_trace);
    end
end

%% Now, for each cell, find if there is a location that it is the only cell been stimulated 

%% 
figure(1)
plot(stim_unique,prob_trace)
%% Now calculate the first spike probability for each cell 

cells_probabilities = zeros([size(trials_locations,1) n_cell]);
for   l = 1:size(trials_locations,1)
    for i_cell = 1:n_cell
    cells_probabilities(l,i_cell)= stim_to_prob(stimuli_size(l,i_cell),stim_unique,prob_trace);
    end
end

%%


%% Cell parameters and background parameters
rng(12242,'twister');
background_rate=1e-4;
v_th_known=15*ones([n_cell,1]);
v_reset_known=-1e4*ones([n_cell,1]);
g_truth = 0.02*ones([n_cell,1]);

funcs.invlink = @invlink_sig;%@(resp) log(1 + exp(resp));%@(x) exp(x);
stim_threshold=10;
time_max =300;

gamma_truth=zeros(n_cell,1);
% gamma_truth = (rand([n_cell 1])<0.5).*(0.5+0.5*rand([n_cell 1]));
gamma_truth(3)=0.8; gamma_truth(10)=0.8;gamma_truth(1)=0.8;gamma_truth(5)=0.8;
gain_truth=0.015+rand([n_cell 1])*0.01;

%% Conduct a few initial random trials targetted at the cells 
n_spots_per_trial = 3;
n_initial_per_plane = round(200*n_cell/n_spots_per_trial/n_planes);
trials_locations =zeros(n_initial_per_plane*n_planes,n_spots_per_trial);
for i_plane=1:n_planes
    for i_trial = 1:n_initial_per_plane
        trials_locations(i_trial + (i_plane-1)*n_initial_per_plane,:)= ...
            randsample(cell_group_list{i_plane},n_spots_per_trial);
    end
end
trials_powers=power_level(3)*ones(n_initial_per_plane*n_planes,1);
%% Need a function to take in the trials locations and others and spit out data
[mpp_temp] = draw_samples(...
    trials_locations, trials_powers, pi_target, background_rate,...
    v_th_known, v_reset_known, g_truth, gain_truth,gamma_truth,...
    current_template,  funcs,    delay_params,stim_threshold,time_max);
%%
mpp=mpp_temp;
%% Variational inference strategy

%histogram(cells_probabilities)
%%
% List the non-identifiable pairs:
% Two cells are non-identifiable if their responses differ only in a few
% trials (<5?)
diff_trials_threshold = 5;

nonidentifiable_pairs=zeros(n_cell,n_cell);
for i_cell = 1:n_cell
    nonid_pairs= find(sum(abs(cells_probabilities-cells_probabilities(:,i_cell)),1)<diff_trials_threshold);
    nonidentifiable_pairs(i_cell,nonid_pairs)=1;
    nonidentifiable_pairs(i_cell,i_cell)=0;
end

%%
disconnected_cell_path = []; %initialize the disconnected cells
n_trial =length(mpp);
eta = 0.001;
%spike_threshold = 40;
S=200;
epsilon=0.01;
changes=1;
n_cell=length(gain_truth);
% psi(x)% digamma function

% Initialize the prior distributions
pi0= 0.7*ones(n_cell,1);alpha0= ones(n_cell,1);beta0 = ones(n_cell,1);

gamma0 = background_rate;

% Initialize the variational family (spike-and-slab with beta distribution)
variational_params=struct([]);
for i_cell = 1:n_cell
    variational_params(i_cell).pi = 0.2;
    variational_params(i_cell).alpha = 1;
    variational_params(i_cell).beta = 1;
end

%---------------------------------------------------%
% Use gradient descent to find the optimal:

sum_of_logs=zeros(S,1);
dqdpi=zeros(n_cell,S);
dqdalpha=zeros(n_cell,S);
dqdbeta=zeros(n_cell,S);

dELBOdpi=zeros(n_cell,S);
dELBOdalpha=zeros(n_cell,S);
dELBOdbeta=zeros(n_cell,S);

temp_rec_pi=zeros(n_cell,1000);
temp_rec_alpha=zeros(n_cell,1000);
temp_rec_beta=zeros(n_cell,1000);
change_history=zeros(1000,1);

iter = 1;
% pickout the identifiable cells
identifiable_index= sum(nonidentifiable_pairs,2)==0;
%%
while changes > epsilon
    v_pi =[variational_params(:).pi]';
    v_alpha = [variational_params(:).alpha]';
    v_beta = [variational_params(:).beta]';
    
    temp_rec_pi(:,iter)=v_pi;
    temp_rec_alpha(:,iter)=v_alpha;
    temp_rec_beta(:,iter)=v_beta;
    
    
    % v_pi=temp_rec_pi(:,iter);
    % v_alpha=temp_rec_alpha(:,iter);
    % v_beta=temp_rec_beta(:,iter);
    
    iter=iter+1;
    for s= 1:S
        % Draw samples from the variational distribution given the current
        % parameters
        gamma_spike = rand(n_cell,1) > v_pi;
        gamma_slab = betarnd(v_alpha,v_beta,[n_cell 1]);
        gamma_sample = gamma_spike.*gamma_slab;
        %bound gamma from 1 to avoid singularity
        gamma_sample=min(gamma_sample, 0.999);
        
        
        % Calculate the loglikelihood given the current samples of gamma
        
        % Need to change how we estimate the likelihood
        
        loglklh_vec = zeros(n_trial,1);
        for i_trial = 1:n_trial
            n_events=length(mpp(i_trial).times);
            [lklh] = calculate_likelihood_sum_bernoulli(n_events,...
                [gamma0; gamma_sample],[1 cells_probabilities(i_trial,:)]');
                loglklh_vec(i_trial)=log(lklh);
        end
        
        loglklh=sum(loglklh_vec);
        %loglklh
        if isinf(loglklh)
            dELBOdpi(:,s)=0;
            dELBOdalpha(:,s)= 0;
            dELBOdbeta(:,s)= 0;
            
        else
            % Calculate the prior probability given the current sample of gamma
            logprior=sum(log(pi0).*(gamma_sample==0)+log(1-pi0).*(gamma_sample>0)+ ...
                (gamma_sample>0).*log(betapdf(gamma_sample,alpha0,beta0)));
            
            % Calculate the probability of the variational distribution given the
            % current sample of gamma
            logvariational=sum(log(max(0.001,v_pi)).*(gamma_sample==0)+...
                log(max(0.001, 1-v_pi)).*(gamma_sample>0)+ ...
                (gamma_sample>0).*log( min(1000,max(0.0001,betapdf(gamma_sample,v_alpha,v_beta)))));
            
            sum_of_logs(s) = loglklh+logprior-logvariational;
            
            % Calculate the gradients of the variational distribution w.r.t. the
            % prior parameters
            dqdpi(:,s)= (gamma_sample==0)./max(0.001,v_pi)-...
                (gamma_sample>0)./max(0.001, 1-v_pi);
            alphaplusbeta = v_alpha+v_beta;
            dqdalpha(:,s)= (gamma_sample>0).*(log(max(0.001,gamma_sample)) + psi(alphaplusbeta)-psi(v_alpha));
            dqdalpha(gamma_sample==0,s)=0;
            dqdbeta(:,s)= (gamma_sample>0).*(log(max(0.001,1-gamma_sample)) + psi(alphaplusbeta)-psi(v_beta));
            dqdbeta(gamma_sample==0,s)=0;
            
            dELBOdpi(:,s)=sum_of_logs(s)*dqdpi(:,s);
            dELBOdalpha(:,s)= sum_of_logs(s)*dqdalpha(:,s);
            dELBOdbeta(:,s)= sum_of_logs(s)*dqdbeta(:,s);
        end
        %dELBOdbeta(:,s)
    end
    % Update the variational parameters using gradient descents
    v_pi_old = v_pi;v_alpha_old = v_alpha;v_beta_old = v_beta;
    v_pi = v_pi+eta*mean(dELBOdpi,2);
    v_alpha = v_alpha+eta*mean(dELBOdalpha,2);
    v_beta = v_beta+eta*mean(dELBOdbeta,2);
    
    v_pi=max(0,min(1,v_pi));
    v_alpha=max(0.001,v_alpha);
    v_beta=max(0.001,v_beta);
    
    
    % Calculate the stopping criteriar
    %changes=sqrt(sum(mean(dELBOdpi,2).^2)+sum(mean(dELBOdalpha,2).^2)+sum(mean(dELBOdbeta,2).^2));
    mean_gamma= (1-v_pi).*v_alpha./(v_alpha+v_beta);
    mean_gamma_old= (1-v_pi_old).*v_alpha_old./(v_alpha_old+v_beta_old);
    
    changes = sqrt(sum((mean_gamma(identifiable_index)-mean_gamma_old(identifiable_index)).^2)/sum(mean_gamma_old(identifiable_index).^2 ));
    changes2=  sum(abs([variational_params(:).pi]'-v_pi))+...
        sum(abs([variational_params(:).alpha]'-v_alpha))+...
        sum(abs([variational_params(:).beta]'-v_beta));
    
    for i_cell = 1:n_cell
        variational_params(i_cell).pi=v_pi(i_cell);
        variational_params(i_cell).alpha=v_alpha(i_cell);
        variational_params(i_cell).beta = v_beta(i_cell);
    end
    change_history(iter)=changes;
    
    fprintf('Iteration: %d; change: %d and %d\n',iter,changes,changes2)
end

% The posterior mean is
mean_gamma= (1-v_pi).*v_alpha./(v_alpha+v_beta)
gamma_truth


%%
% Filter the disconneced neurons (very low gamma mean with small variance)

variance_gamma = (1-v_pi).*(mean_gamma.^2+ v_alpha.*v_beta./(v_alpha+v_beta).^2./(v_alpha+v_beta+1));

disconnected_cells = (mean_gamma<0.1) & (variance_gamma <0.01);
disconnected_cells = disconnected_cells | disconnected_cells;
% how to incorporate previous info on disconnected cells?
% add the list of disconnected cell to the previous one?
disconnected_cell_path = [disconnected_cell_path disconnected_cells];


%% Calculate the posterior variances
post_variances=variance_gamma; % useful for calculating the score for each locations later
post_variances(sum(nonidentifiable_pairs,1)>0)=1000;

% set the variances for the unidentifiable pairs to be infinity


%% Now predict stim reactions
% Evaluate the expected output on the target locations

% Set these parameters as the prior parameters:
pi0=[variational_params(:).pi]';
alpha0=[variational_params(:).alpha]';
beta0=[variational_params(:).beta]';


epsilon=0.01;
params_by_target=struct([]);
score_target = zeros(size(target_locations,1),1);
for i_target= 1:size(target_locations,1)
    % calculate the size of stimulations, and find the cells stimulated
    stimuli_temp = pi_target(:,i_target).*power_level(3);
    temp_idx= stimuli_temp>spike_threshold; %cells_stimulated_temp
    
    % If the selected cells are identifiable in the previous data, we will
    % consider this locations:
    if sum(nonidentifiable_pairs(min(find(temp_idx)),find(temp_idx)))>0
        % meaning that the new stimulated cells are non-identifiable
        score_target(i_target)=0;
        
    else
        
        % Narrow our discussion to cells that are stimulated in this trials
        pi0_local=[variational_params(temp_idx).pi]';
        % bound pi0_local away from singularity:
        pi0_local = min(0.9,max(0.1,pi0_local));
        
        alpha0_local=[variational_params(temp_idx).alpha]';
        beta0_local=[variational_params(temp_idx).beta]';
        
        
        
        % Variational updates:
        
        
        
        temp_rec_pi=zeros(sum(temp_idx),1000);
        temp_rec_alpha=zeros(sum(temp_idx),1000);
        temp_rec_beta=zeros(sum(temp_idx),1000);
        change_history=zeros(1000,1);
        
        changes=1;
        v_pi_old =pi0_local;
        v_alpha_old = alpha0_local;
        v_beta_old = beta0_local;
        iter=1;
        while changes > epsilon
            v_pi =v_pi_old;
            v_alpha = v_alpha_old;
            v_beta = v_beta_old;
            
            temp_rec_pi(:,iter)=v_pi;
            temp_rec_alpha(:,iter)=v_alpha;
            temp_rec_beta(:,iter)=v_beta;
            
            % v_pi=temp_rec_pi(:,iter);
            % v_alpha=temp_rec_alpha(:,iter);
            % v_beta=temp_rec_beta(:,iter);
            iter=iter+1;
            dELBOdpi_temp=zeros(sum(temp_idx),S);
            dELBOdalpha_temp=zeros(sum(temp_idx),S);
            dELBOdbeta_temp=zeros(sum(temp_idx),S);
            gsample=zeros(1,S);
            for s= 1:S
                gamma_spike = rand(sum(temp_idx),1) > v_pi;
                gamma_slab = betarnd(v_alpha,v_beta,[sum(temp_idx) 1]);
                gamma_sample = gamma_spike.*gamma_slab;
                %bound gamma from 1 to avoid singularity
                gamma_sample=min(gamma_sample, 0.999);
                gsample(s)=gamma_sample;
                
                % Calculate the gradient given the the prob_no_event, and prob_event
                prob_no_event= (1-gamma0)*prod(1-gamma_sample);
                prob_event=1-prob_no_event;
                
                % Calculate the prior probability given the current sample of gamma
                logprior=sum(log(pi0_local).*(gamma_sample==0)+log(1-pi0_local).*(gamma_sample>0)+ ...
                    (gamma_sample>0).*log( min(1000,max(0.0001,betapdf(gamma_sample,v_alpha,v_beta)))));
                
                % Calculate the probability of the variational distribution given the
                % current sample of gamma
                logvariational=sum(log(max(0.001,v_pi)).*(gamma_sample==0)+...
                    log(max(0.001, 1-v_pi)).*(gamma_sample>0)+ ...
                    (gamma_sample>0).*log( min(1000,max(0.0001,betapdf(gamma_sample,v_alpha,v_beta)))));
                
                % Calculate the gradients of the variational distribution w.r.t. the
                % prior parameters
                dqdpi_temp= (gamma_sample==0)./max(0.001,v_pi)-...
                    (gamma_sample>0)./max(0.001, 1-v_pi);
                alphaplusbeta = v_alpha+v_beta;
                dqdalpha_temp= (gamma_sample>0).*(log(max(0.001,gamma_sample)) + ...
                    psi(alphaplusbeta)-psi(v_alpha));
                dqdalpha_temp(gamma_sample==0)=0;
                dqdbeta_temp= (gamma_sample>0).*(log(max(0.001,1-gamma_sample)) + ...
                    psi(alphaplusbeta)-psi(v_beta));
                dqdbeta_temp(gamma_sample==0)=0;
                
                
                % Calculate the loglikelihood given the current samples of gamma
                if prob_no_event>0 % avoid singularity
                    loglklh_0=log(prob_no_event);
                    sum_of_logs_0 = loglklh_0+logprior-logvariational;
                    dELBOdpi_temp(:,s)=sum_of_logs_0*dqdpi_temp*prob_no_event;
                    dELBOdalpha_temp(:,s)=sum_of_logs_0*dqdalpha_temp*prob_no_event;
                    dELBOdbeta_temp(:,s)=sum_of_logs_0*dqdbeta_temp*prob_no_event;
                end
                
                if prob_event > 0 % avoid singularity
                    loglklh_1=log(prob_event);
                    sum_of_logs_1 = loglklh_1+logprior-logvariational;
                    dELBOdpi_temp(:,s)=dELBOdpi_temp(:,s)+sum_of_logs_0*dqdpi_temp*prob_event;
                    dELBOdalpha_temp(:,s)=dELBOdalpha_temp(:,s)+sum_of_logs_0*dqdalpha_temp*prob_event;
                    dELBOdbeta_temp(:,s)=dELBOdbeta_temp(:,s)+sum_of_logs_0*dqdbeta_temp*prob_event;
                end
                
            end
            
            v_pi = v_pi+eta*mean(dELBOdpi_temp,2);
            v_alpha = v_alpha+eta*mean(dELBOdalpha_temp,2);
            v_beta = v_beta+eta*mean(dELBOdbeta_temp,2);
            
            v_pi=max(0,min(1,v_pi));
            v_alpha=max(0.001,v_alpha);
            v_beta=max(0.001,v_beta);
            
            % Calculate the stopping criteriar
            mean_gamma= (1-v_pi).*v_alpha./(v_alpha+v_beta);
            mean_gamma_old= (1-v_pi_old).*v_alpha_old./(v_alpha_old+v_beta_old);
            
            changes = sqrt(sum((mean_gamma-mean_gamma_old).^2)/sum(mean_gamma_old.^2 ));
            change_history(iter)=changes;
            v_pi_old = v_pi;
            v_alpha_old = v_alpha;
            v_beta_old = v_beta;
            
            fprintf('Iteration: %d; change: %d\n',iter,changes)
        end
        
        % Calculate the score of this target as the changes in the posterior
        % variances
        mean_gamma= (1-v_pi).*v_alpha./(v_alpha+v_beta);
        new_variance=(1-v_pi).*(mean_gamma.^2+ v_alpha.*v_beta./(v_alpha+v_beta).^2./(v_alpha+v_beta+1));
        
        score_target(i_target) = sum(abs(post_variances(temp_idx)-new_variance));
        
    end
    
end
%% Now select the trials based on the scores


%% Testing code:


% Find out how many cell are stimulated in each trial
spike_threshold = 40;
n_cell=length(gain_truth);
stimuli_size=zeros(size(trials_locations,1),n_cell);
for l = 1:size(trials_locations,1)
    for m = 1:size(trials_locations,2)
        if isnan(trials_locations(l,m))
        else
            stimuli_size(l,:) = stimuli_size(l,:)+...
                (pi_target(:,trials_locations(l,m)).*trials_powers(l))';
        end
    end
end

cells_stimulated = zeros([size(trials_locations,1) n_cell]);
for   l = 1:size(trials_locations,1)
    cells_stimulated(l,: )= stimuli_size(l,:)>spike_threshold;
end

% need a vector to represent probability of events from other sources
% the probability depends on the number of other stimulated cells, and the
% distribution of gamma
% can be an estimates from simulation
sparsity_level = 0.05;
prob_spikes_other = zeros(21,1); % 0 to 20 cells
for i_events=1:length(prob_spikes_other)
    prob_no_spikes = (1-time_max*background_rate)*(1-sparsity_level)^(i_events-1);
    prob_spikes_other(i_events) = 1-prob_no_spikes;
end

%%
% We can calculate the likelihood across different values of gamma
% ...find MLE of gamma, and check its variance?

gamma_low = 0.2;
gamma_up = 0.8;
prob_gamma_greater_than_lb = zeros(n_cell,1);
prob_gamma_smaller_than_ub = zeros(n_cell,1);

for i_cell = 1:n_cell
    related_trials = find(cells_stimulated(:,i_cell));
    prob_low=1;prob_up=1;
    for i_trial = 1:length(related_trials)
        num_stimulated_cells=sum(cells_stimulated(related_trials(i_trial),:));
        prob_low_temp = prob_spikes_other(num_stimulated_cells+1)*(1-gamma_low) +gamma_low;
        prob_up_temp = prob_spikes_other(num_stimulated_cells+1)*(1-gamma_up) +gamma_up;
        if isempty(mpp_temp(related_trials(i_trial)).times)
            prob_low = prob_low*(1-prob_low_temp);
            prob_up = prob_up*(1-prob_up_temp);
        else
            prob_low = prob_low*(prob_low_temp);
            prob_up = prob_up*(prob_up_temp);
        end
    end
    prob_gamma_greater_than_lb(i_cell)=prob_low;
    prob_gamma_smaller_than_ub(i_cell)=prob_up;
end
% We can simulate the null distribution given the number of cells
% stimulated in each trial

