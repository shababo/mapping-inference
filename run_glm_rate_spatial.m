function glm_out = run_glm_rate_spatial(results_grid)

num_locations = numel(results_grid);

num_trials = num_locations * length(results_grid{1,1});

% stim_mat = zeros(num_trials,num_locations);
stim_mat = zeros(num_trials,num_locations + 1);
stim_mat(:,1) = 1;
response = zeros(num_trials,1);
response_map = zeros(size(results_grid));

trial_i = 1;

for i = 1:size(results_grid,1)
    for j = 1:size(results_grid,2)
        
        for k = 1:length(results_grid{i,j})
            
            stim_mat(trial_i,(i-1)*size(results_grid,2) + j + 1) = 1;
            psth = histcounts(results_grid{i,j}(k).times,0:20:800);
            psth = smoothts(psth,'g',10,2);
            response(trial_i) = sum(psth);
            response_map(i,j) = response_map(i,j) + sum(psth);
            trial_i = trial_i + 1;
        end
    end
end

% [~, ~, glm_out] = glmfit(stim_mat,response,'poisson');
glm_out = cvglmnet(stim_mat,response,'gaussian');
% glm_out = [];
figure; imagesc(response_map)