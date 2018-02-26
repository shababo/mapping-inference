function [loglklh] = update_likelihood_dev2(...
    gamma_sample_mat,gain_sample_mat,stim_size,mpp,...
    prob_trace_full,minimum_stim_threshold,stim_scale,background_rate,...
    relevant_trials,lklh_func,spike_curves)

n_cell=size(gamma_sample_mat,1);
S=size(gamma_sample_mat,2);
n_trial=length(mpp);
% n_grid=size(prob_trace_full,2);
n_grid=300;
t_grid= 1:n_grid;
loglklh=zeros(n_cell,S);
for s=1:S
    gamma_sample=gamma_sample_mat(:,s);gain_sample=gain_sample_mat(:,s);
    loglklh_vec = zeros(n_trial,1);
    for  i_trial = 1:n_trial
        
        stim_temp =stim_size(i_trial,:)';
        effective_stim= stim_temp.*gain_sample;
        stimulated_cells = find(effective_stim>minimum_stim_threshold);
        effective_stim=effective_stim(stimulated_cells );
        
        % find the index of the actual stimulation in the vector:
        if ~isempty(stimulated_cells)
            stim_index=zeros(length(stimulated_cells),1);
            for i_stim = 1:length(stimulated_cells)
            [~, stim_index(i_stim)]=min(abs(effective_stim(i_stim) - spike_curves.current));
            end
        end
        
        %stim_index=min(n_grid,max(1,round(effective_stim*stim_scale)));
%       prob_this_trial= (gamma_sample(stimulated_cells)*ones(1,n_grid)).*prob_trace_full(stim_index,:);
%       prob_this_trial=[background_rate*ones(1, size(prob_this_trial,2)); prob_this_trial];
        if ~isempty(stimulated_cells)
            prob_this_trial=zeros(length(stimulated_cells),n_grid);
            for i_stim = 1:length(stimulated_cells)
                prob_this_trial(i_stim,:)=...
                    normpdf(t_grid,spike_curves.mean(stim_index(i_stim)),spike_curves.sd(stim_index(i_stim)));
            end
            prob_this_trial=[background_rate*ones(1, n_grid); prob_this_trial];
        else
            prob_this_trial=[background_rate*ones(1, n_grid)];
        end
        
        [loglklh_vec(i_trial)]=  lklh_func(mpp(i_trial),prob_this_trial);
       
    end
    for i_cell = 1:n_cell
        loglklh(i_cell,s)=sum(loglklh_vec(relevant_trials{i_cell}));
    end
end

end
