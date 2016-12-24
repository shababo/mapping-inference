function [glm_out, response_map] = run_glm_num_events_spatial_mapest(map_est,time_window,resp_only)

num_locations = numel(map_est);

num_trials = num_locations * length(map_est{1,1});

% stim_mat = zeros(num_trials,num_locations);
stim_mat = zeros(num_trials,num_locations + 1);
stim_mat(:,1) = 1;
response = zeros(num_trials,1);
response_map = zeros(size(map_est));
trial_i = 1;

for i = 1:size(map_est,1)
    for j = 1:size(map_est,2)
        these_map_ests = map_est{i,j};
        for k = 1:length(map_est{i,j})
            
            stim_mat(trial_i,(i-1)*size(map_est,2) + j + 1) = 1;
%             psth = histcounts(results_grid{i,j}(k).times,1:1:2000);
%             psth = smoothts(psth,'g',100,20);
%             [pks, locs] = findpeaks(psth(150:800),...
%                 'MinPeakHeight',.25*std(psth),'MinPeakDistance',20);
%             [~,map_ind] = min(results_grid{i,j}(k).obj);
%             map_sample = ...
%                 truncate_samples(results_grid{i,j}(k),[map_ind map_ind]);
            this_map_est = these_map_ests{k};
            [this_map_est.times, sort_inds] = sort(this_map_est.times);
            this_map_est.amp = this_map_est.amp(sort_inds);
            this_map_est.amp(find(diff(this_map_est.times) < 60) + 1) = [];
            this_map_est.times(find(diff(this_map_est.times) < 60) + 1) = [];
            
            response(trial_i) = sum(this_map_est.times >= time_window(1) & this_map_est.times <= time_window(2) & this_map_est.amp >= 0);
%             if response(trial_i) > 5
%                 response(trial_i) = 5;
%             end
            response_map(i,j) = response_map(i,j) + response(trial_i);
            trial_i = trial_i + 1;
            
            
        end
%         response_map(i,j) = response_map(i,j)/length(map_est{i,j});
    end
end

if ~resp_only
    % [~, ~, glm_out] = glmfit(stim_mat,response,'poisson');
    options.alpha = 1.0;
    options.nfolds = 3;
    glm_out = cvglmnet(stim_mat,response,'poisson',options);
else
    glm_out = 0;
end
% figure;
% imagesc(response_map)