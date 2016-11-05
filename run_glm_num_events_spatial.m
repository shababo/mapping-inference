function glm_out = run_glm_num_events_spatial(num_events_grid)

num_locations = numel(num_events_grid);

num_trials = num_locations * length(num_events_grid{1,1});

% stim_mat = zeros(num_trials,num_locations);
stim_mat = zeros(num_trials,num_locations + 1);
stim_mat(:,1) = 1;
response = zeros(num_trials,1);

trial_i = 1;

for i = 1:size(num_events_grid,1)
    for j = 1:size(num_events_grid,2)
        
        for k = 1:length(num_events_grid{i,j})
            
            stim_mat(trial_i,(i-1)*size(num_events_grid,2) + j + 1) = 1;
            response(trial_i) = num_events_grid{i,j}(k);
            
            trial_i = trial_i + 1;
        end
    end
end

% [~, ~, glm_out] = glmfit(stim_mat,response,'poisson');
options.alpha = .5;
glm_out = cvglmnet(stim_mat,response,'poisson',options);