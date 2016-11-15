function glm_out = run_glm_num_events_spatial(results_grid)

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
%             psth = histcounts(results_grid{i,j}(k).times,1:1:2000);
%             psth = smoothts(psth,'g',100,20);
%             [pks, locs] = findpeaks(psth(150:800),...
%                 'MinPeakHeight',.25*std(psth),'MinPeakDistance',20);
            [~,map_ind] = min(results_grid{i,j}(k).obj);
            map_sample = ...
                truncate_samples(results_grid{i,j}(k),[map_ind map_ind]);
            response(trial_i) = length(map_sample.times);
            if response(trial_i) > 5
                response(trial_i) = 5;
            end
            response_map(i,j) = response_map(i,j) + response(trial_i);
            trial_i = trial_i + 1;
            
            
        end
        response_map(i,j) = response_map(i,j)/length(results_grid{i,j});
    end
end

% [~, ~, glm_out] = glmfit(stim_mat,response,'poisson');
% options.alpha = 1.0;
glm_out = cvglmnet(stim_mat,response,'poisson');
figure;
imagesc(response_map)