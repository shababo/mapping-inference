%gamma_ini = 0.5;
gamma_current = overall_connectivity;

mu_current = overall_mark;
%sigma_current = evoked_params.sigma_a*ones(n_cell,1);
sigma_current = 1*ones(n_cell_local,1);

mu_m_hyper =overall_mark;
mu_v_hyper = ones(n_cell_local,1); % inverse of prior variance 
sigma_alpha = ones(n_cell_local,1);
sigma_beta = ones(n_cell_local,1);

gamma_alpha = overall_connectivity;
gamma_beta = ones(n_cell_local,1);

f_background = bg_params.firing_rate/1000; 
w_background = 1/(bg_params.a_max -bg_params.a_min); % size distribution for background events
%% 
if exact_crossing == 1 
    estimated_intensity = E_intensity;
else 
    estimated_intensity = M_intensity;
end
%% Precalculation:
expected_all = zeros(i_trial,i_cell);
for i_trial = 1:n_trial
    for i_cell = 1:n_cell_local
        expected_all(i_trial,i_cell)=sum(estimated_intensity{i_trial, i_cell} );
    end
end
%% 
mini_factor = n_trial/n_trial_update;
        
%% Run Gibbs sampler
gamma_samples = zeros(n_gibbs_sample, n_cell_local);
mu_samples = zeros(n_gibbs_sample, n_cell_local);
sigma_samples = zeros(n_gibbs_sample, n_cell_local);
soft_assignments_samples = cell(n_gibbs_sample,1);
soft_assignments_current = cell(n_trial_update,1);
i_sample = 0;
i_counter = 0;
while i_sample < n_gibbs_sample
    % Gibbs sampler
    %events_precell = cell(n_cell,1);
    events_precell = cell(n_cell_local,1);
    for i_cell = 1:n_cell_local
        events_precell{i_cell} = [];
    end
    
    % Obtain the list of stimulated cells 
    evoked_cell_batch = [];
    chosen_trials_index = randsample(n_trial, n_trial_update);
    for i_trial_index = 1:n_trial_update
        i_trial = chosen_trials_index(i_trial_index);
        evoked_cell_batch = [evoked_cell_batch evoked_cell{i_trial}];
    end
    evoked_cell_batch = unique(evoked_cell_batch);
    evoked_cell_batch = evoked_cell_batch(2:end); % Dropping the spontaneous events
    
    n_cell_evoked = length(evoked_cell_batch);
    
    % Draw assignments
    % Pick a few trials
    for i_trial_index = 1:n_trial_update
        i_trial = chosen_trials_index(i_trial_index);
        n_events = length(mpp_new(i_trial).event_times);
        evoked_cell_index = evoked_cell{i_trial};
        firing_rates = zeros(length(evoked_cell_index),1);
        size_rates = zeros(length(evoked_cell_index),1);
        soft_assignments_current{i_trial_index}  = ones(n_events, length(evoked_cell_index))/length(evoked_cell_index);
        
        for i_event = 1:n_events
            this_event_time = mpp_new(i_trial).event_times(i_event);
            this_event_size = mpp_new(i_trial).amplitudes(i_event);
            [~, i_t]=min(abs(t_vect - this_event_time));
                
            for i_index = 1:length(evoked_cell_index)
                 i_cell=evoked_cell_index(i_index);
                if i_cell== 0
                    firing_rates(i_index)  = f_background;
                    size_rates(i_index) = w_background*(this_event_size>bg_params.a_min)*...
                        (this_event_size<bg_params.a_max); 
                else 
                     firing_rates(i_index)  = gamma_current(i_cell)*estimated_intensity{i_trial,i_cell}(i_t);
                    size_rates(i_index) = normpdf(this_event_size,mu_current(i_cell),sigma_current(i_cell));
                end
            end
            
            % Draw assignments given the estimated rates
            %chances = firing_rates.*size_rates;
            chances = firing_rates;
            if sum(chances)==0
               chances(1) = 1; 
            end
            chances = chances/sum(chances);
            %fprintf('%d', sum(chances));
            soft_assignments_current{i_trial_index}(i_event,:) = chances; 
        end
        for i_index = 2:length(evoked_cell_index)
                i_cell = evoked_cell_index(i_index);
                events_precell{i_cell} = [events_precell{i_cell} [mpp_new(i_trial).amplitudes; soft_assignments_current{i_trial_index}(:,i_index)']];
        end 
    end
    
    %-------------------------------------------------------%
    % Draw mu and sigma 
    if sigma_unknown == 1
        %Not available...
    else
        for i = 1:n_cell_evoked
            i_cell= evoked_cell_batch(i);
            % Sum of the weights:
            if  sum(size(events_precell{i_cell}))> 1
            weighted_sum = sum(events_precell{i_cell},2);
            weighted_sum(1) = sum(events_precell{i_cell}(1,:).*events_precell{i_cell}(2,:));
            if weighted_sum(2) > 0
                post_mean = (mu_m_hyper(i_cell)*mu_v_hyper(i_cell)+mini_factor*weighted_sum(1))/...
                    (mu_v_hyper(i_cell)+ mini_factor*weighted_sum(2));
                post_v = mu_v_hyper(i_cell)+mini_factor*weighted_sum(2);
            else
                post_mean = mu_m_hyper(i_cell);
                post_v = mu_v_hyper(i_cell);
            end
            % Draw mu:
            mu_current(i_cell) = normrnd(post_mean,sigma_current(i_cell)/sqrt(post_v) );
            else
                % Do nothing
            end
        end
    end
    
    % Draw gamma
    % We take an approximated approach:
    %   1) Calculate the expected number of events
    %   2) Update the beta and gamma given the number of success
    %   3) Draw new gammas
    for i = 1:n_cell_evoked
            i_cell= evoked_cell_batch(i);
            if sum(size(events_precell{i_cell}))> 1
            weighted_sum = sum(events_precell{i_cell},2);
            weighted_sum(1) = sum(events_precell{i_cell}(1,:).*events_precell{i_cell}(2,:));
            expected_events = sum( expected_all(chosen_trials_index,i_cell));
            post_alpha = gamma_alpha(i_cell)+mini_factor*weighted_sum(2);
            post_beta = gamma_beta(i_cell) + mini_factor*(expected_events-weighted_sum(2));
            % Draw gamma:
            if expected_events < weighted_sum(2)
                gamma_current(i_cell)=1;
            else
                gamma_current(i_cell) = betarnd(post_alpha,post_beta);
            end
            if isnan( gamma_current(i_cell))
                fprintf('Wrong')
                fprintf('%d\n', i_cell)
                break
            end
        else
            expected_events = sum( expected_all(chosen_trials_index,i_cell));
            post_alpha = gamma_alpha(i_cell)+0;
            post_beta = gamma_beta(i_cell) + mini_factor*expected_events-0;
            % Draw gamma:
            if expected_events < weighted_sum(2)
                gamma_current(i_cell)=1;
            else
                gamma_current(i_cell) = betarnd(post_alpha,post_beta);
            end
            end        
    end
    
    %----------------------------------------------%
    % Take independent samples
    i_counter = i_counter +1;
    if i_counter > n_burnin 
        if mod(i_counter-1,n_skip) == 0  
            i_sample = i_sample +1;
            gamma_samples(i_sample,:) = gamma_current;
            mu_samples(i_sample,:) = mu_current;
            sigma_samples(i_sample,:) = sigma_current;
            
            soft_assignments_samples{i_sample}  = soft_assignments_current;  
        fprintf('%d\n',i_sample);
        end
    else
        fprintf('burnin %d\n',i_counter);
    end
    
end





