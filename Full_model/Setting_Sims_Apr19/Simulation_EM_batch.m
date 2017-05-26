%% Run the EM algorithm
gamma_samples = gamma_old;
mu_samples = mu_old;
sigma_samples = sigma_old;

soft_assignments_samples = cell(1,1);
soft_assignments = cell(n_trial_update,1);



chosen_trials_index = 1:n_trial_temp;

evoked_cell_batch = unique([evoked_cell{chosen_trials_index}]);
evoked_cell_batch = evoked_cell_batch(2:end); % Drop the spontaneous events
n_cell_evoked = length(evoked_cell_batch);

count_trials = zeros(n_cell_evoked,1);
count_events = zeros(n_cell_evoked,1);
for i_trial_index = 1:n_trial_update
    i_trial = chosen_trials_index(i_trial_index);
    n_events= length( mpp_new(i_trial).amplitudes);
    if length(evoked_cell{i_trial})>1
    count_trials(evoked_cell{i_trial}(2:end)) = count_trials(evoked_cell{i_trial}(2:end))+1;
    count_events(evoked_cell{i_trial}(2:end))  = count_events(evoked_cell{i_trial}(2:end)) +...
        n_events;
    end
    soft_assignments{i_trial_index}  = ones(n_events, length(evoked_cell{i_trial}))/length(evoked_cell{i_trial});
    
end
events_precell = cell(n_cell_evoked,1);
trials_precell = cell(n_cell_evoked,1);
times_precell = cell(n_cell_evoked,1);

for i_ind = 1:n_cell_evoked
    i_cell = evoked_cell_batch(i_ind);
    events_precell{i_cell} = zeros(2,count_events(i_cell));
    trials_precell{i_cell} = zeros(count_trials(i_cell),1);
    times_precell{i_cell} = zeros(2,count_events(i_cell));
end

count_trials(:)=0;
count_events(:)=0;
for i_trial_index = 1:n_trial_update
    i_trial = chosen_trials_index(i_trial_index);
    evoked_cell_index = evoked_cell{i_trial};
    
    this_trial_time = mpp_new(i_trial).event_times;
    this_trial_amps = mpp_new(i_trial).amplitudes;
    n_events=length(this_trial_time);
    if length(evoked_cell_index)>1
        for i_index = 2:length(evoked_cell_index)
            i_cell = evoked_cell_index(i_index);
            events_precell{i_cell}(1,count_events(i_cell) + (1:n_events)) = this_trial_amps;
            %events_precell{i_cell}(2,count_events(i_cell) + (1:n_events)) = soft_assignments{i_trial_index}(:,i_index);
            times_precell{i_cell}(1,count_events(i_cell) + (1:n_events))  =this_trial_time;
            times_precell{i_cell}(2,count_events(i_cell) + (1:n_events))  =i_trial;
            trials_precell{i_cell}(count_trials(i_cell)+1) = i_trial;
            count_trials(i_cell) = count_trials(i_cell) +1;
            count_events(i_cell) =count_events(i_cell)+ n_events;
        end
    end
end



normalized_change = convergence_epsilon + 1;
num_iter = 0;

while (normalized_change > convergence_epsilon) & (num_iter < maxit)
    
    num_iter = num_iter+1;
    % Draw assignments
    % Pick a few trials
    temp_ind = 0;
    count_events(:)=0;
    
    for i_trial_index = 1:n_trial_update
        i_trial = chosen_trials_index(i_trial_index);
        n_events = length(mpp_new(i_trial).event_times);
        evoked_cell_index = evoked_cell{i_trial};
        firing_rates = zeros(length(evoked_cell_index),1);
        size_rates = zeros(length(evoked_cell_index),1);
       
          this_trial_time = mpp_new(i_trial).event_times;
          this_trial_amps = mpp_new(i_trial).amplitudes;
         for i_event = 1:n_events
            this_event_time =this_trial_time(i_event);
            this_event_size =this_trial_amps(i_event);
            [~, i_t]=min(abs(t_vect - this_event_time));
            for i_index = 1:length(evoked_cell_index)
                i_cell=evoked_cell_index(i_index);
                if i_cell== 0
                    firing_rates(i_index)  = f_background;
                    size_rates(i_index) = ...
                        exp(-(log(this_event_size)-mean_background).^2/(2*sigma_background^2))/(sqrt(2*pi)*sigma_background);
                else
                    firing_rates(i_index)  = gamma_old(i_cell)*estimated_intensity{i_trial,i_cell}(i_t);
                    size_rates(i_index) = ...
                        exp(-(this_event_size-mu_old(i_cell)).^2/(2*sigma_old(i_cell)^2))/(sqrt(2*pi)*sigma_old(i_cell));
                end
            end
            % Draw assignments given the estimated rates
            chances = firing_rates.*size_rates;
            if sum(chances)==0
               chances(1) = 1; 
            end
            chances = chances/sum(chances);
            %fprintf('%d', sum(chances));
            soft_assignments{i_trial_index}(i_event,:) = chances; 
            if isnan( chances(1))
                fprintf('Wrong')
                fprintf('%d\n', i_cell)
                temp_ind = 1;
                break
            end
           
        end
               if length(evoked_cell_index)>1
        for i_index = 2:length(evoked_cell_index)
                i_cell = evoked_cell_index(i_index);
                events_precell{i_cell}(2,count_events(i_cell) + (1:n_events)) = soft_assignments{i_trial_index}(:,i_index);
                count_events(i_cell) =count_events(i_cell)+ n_events;
        end 
               end
        
    end
        if temp_ind == 1
            break
        end
    
    % The M-step:
    gamma_current = gamma_old;
    mu_current = mu_old;
    sigma_current = sigma_old;
   
    %-------------------------------------------------------%
    % Update mu and sigma
    for i = 1:n_cell_evoked
        i_cell= evoked_cell_batch(i);
        % Sum of the weights:
        if  sum(size(events_precell{i_cell}))> 1
            weighted_sum = sum(events_precell{i_cell},2);
            if weighted_sum(2) > 1
                weighted_sum(1) = sum(events_precell{i_cell}(1,:).*events_precell{i_cell}(2,:));
                mu_current(i_cell) = weighted_sum(1)/weighted_sum(2);
                if sigma_unknown == 1
                    weighted_sigma = sum(events_precell{i_cell}(2,:).*((events_precell{i_cell}(1,:)-mu_current(i_cell)).^2));
                    sigma_current(i_cell) = sqrt(weighted_sigma/weighted_sum(2));
                else
                    % Do nothing since sigma is fixed
                end
            else
                % Do nothing...not enough assigned events
            end
        else
            % Do nothing...not enough times stimulated
        end
    end
    % Update gamma
    for i = 1:n_cell_evoked
        i_cell= evoked_cell_batch(i);
        if sum(size(events_precell{i_cell}))> 1
            weighted_sum = sum(events_precell{i_cell},2);
            expected_events = sum( expected_all(trials_precell{i_cell},i_cell));
            gamma_current(i_cell) = min(1,weighted_sum(2)/expected_events);
            % gamma shouldn't be larger than one
            if isnan( gamma_current(i_cell))
                fprintf('Wrong')
                fprintf('%d\n', i_cell)
                break
            end
        else
            % Do nothing since not stimulated enough times
        end
    end
    % Evaluate the convergence
    normalized_change = norm(gamma_current - gamma_old)/(norm(gamma_old)+1) + norm(mu_current - mu_old)/(norm(mu_old)+1)+...
        norm(sigma_current - sigma_old)/(norm(sigma_old)+1);
    
    % Update the paremeters
    gamma_old = gamma_current;
    mu_old = mu_current;
    sigma_old = sigma_current;

    %----------------------------------------------%
    % Record the solution path
    if num_iter == 1
        gamma_samples = gamma_current;
        mu_samples = mu_current;
        sigma_samples = sigma_current;
    else
        gamma_samples = [gamma_samples gamma_current];
        mu_samples = [mu_samples mu_current];
        sigma_samples = [sigma_samples sigma_current];
        
    end
    %soft_assignments_samples{num_iter} = soft_assignments;
    
   
    
    fprintf('Iteration: %d, normalized changes %d\n',num_iter, normalized_change);
end

