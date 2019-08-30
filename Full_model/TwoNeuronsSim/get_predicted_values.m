function [new_trials] = get_predicted_values(neurons,  new_trials, prior_info, params)
%% Outline
% - Draw samples of all relevant parameters from their posterior
% distributions
% - Draw samples of the new locations from the conditional distributions
% (derived based on the prior)
% - Draw spike time and event time given each sample
% - Summarize the posterior samples (parameters and event times)

% params.plot_only
% params.MC_params.sample_size = 50;
n_cell = length(neurons);
Tmax=prior_info.induced_intensity.event_time_max;
time_factor = prior_info.induced_intensity.time_factor;
timepoints = 1:Tmax;
min_dist = prior_info.GP_params.min_dist;
%% Extract posterior information from neurons
clear('posterior_params')
for i_cell = 1:n_cell
    posterior_params(i_cell)=neurons(i_cell).params(end);
end

%% Draw samples from posterior distributions
S= params.MC_params.sample_size;
posterior_samples = cell([S 1]);
for s =1:S
    [posterior_samples{s},~] = draw_samples_from_var_dist(posterior_params);
end
%% If the data contain only single location (no shapes are fitted)
% One can just call the get_fitted_values() function
if ~(isfield(neurons(i_cell).params(end),'shapes') |isfield(neurons(i_cell).params(end),'xy'))
    
    [~,new_trials]=classify_locations(posterior_params,neurons,new_trials,prior_info);
    
    [fitted_trials]=get_fitted_values(neurons, new_trials, prior_info,params);
    for i = 1:length(fitted_trials)
        new_trials(i).predicted=fitted_trials(i).fitted;
    end
else
    
    %% Obtain the conditional distribution
    % The distribution of the shape values at the new locations, conditioning on the stimulated locations
    n_trials = length(new_trials);
    
    % Classify the locations in the new trials
    [new_shape_params,new_trials]=classify_locations(posterior_params,neurons,new_trials,prior_info);
    
    % Calculate the conditional distribution
    [new_shape_params]=get_new_shape_conditional(posterior_params,new_shape_params,prior_info);
    
    
    % Draw samples from the conditional distribution
    %if size(new_shape_params,1)>0
    new_shape_samples = cell([S 1]);
    for s=1:S
        [new_shape_samples{s}]=draw_samples_from_shape_conditional(new_shape_params,posterior_samples{s});
    end
    %end
    
    % else
    %     new_shape_samples=cell([S 1]);
    %     emp_struct=struct;
    %     for i_cell = 1:n_cell
    %         emp_struct(i_cell).shapes=[];
    %     end
    %     for s=1:S
    %         new_shape_samples{s}=emp_struct;%draw_samples_from_shape_conditional(new_shape_params,posterior_samples{s});
    %     end
    % end
    %
    
    %% Calculate the fitted intensities (as posterior averages) on new trials
    bg_tmp =cellfun(@(x) x.background,posterior_samples);
    background_post=mean(bg_tmp);
    
    for i_trial = 1:n_trials
        new_trials(i_trial).predicted=struct;
        new_trials(i_trial).predicted.intensity=struct;
        new_trials(i_trial).predicted.intensity.spike= background_post*ones(1,length(timepoints));
        new_trials(i_trial).predicted.intensity.event=background_post*ones(1,length(timepoints));
        new_trials(i_trial).predicted.PR=1;
        new_trials(i_trial).predicted.source=0;
        new_trials(i_trial).predicted.timepoints=timepoints;
        new_trials(i_trial).predicted.time_factor=time_factor;
        new_trials(i_trial).predicted.stim=0;
    end
    
    for i_cell = 1:n_cell
        
        PR_tmp =cellfun(@(x) x.PR,posterior_samples);
        PR_post=mean(PR_tmp);
        
        for i_trial =1:n_trials
            this_trial=new_trials(i_trial);
            %         rel_pos=neurons(i_cell).params(end).shapes.locations - ...
            %             ones(size(neurons(i_cell).params(end).shapes.locations,1),1)*(this_trial.locations-neurons(i_cell).location);
            %         i_shape=find(sum( (rel_pos.^2)')==0);
            intensity_records=struct;
            intensity_records.spike= zeros(S,length(timepoints));
            intensity_records.event= zeros(S,length(timepoints));
            stim_records=zeros(S,1);
            
            for s=1:S
                this_sample=posterior_samples{s}(i_cell);
                for i_loc = 1:size(new_trials(i_trial).locations,1)
                    cell_and_pos=new_trials(i_trial).cell_and_pos{i_loc};
                    cell_and_status=new_trials(i_trial).cell_and_status{i_loc};
                    
                    stim=0;
                    %                 this_trial.power_levels*this_sample.gain*this_sample.shapes(i_shape);
                    if ~isempty(cell_and_pos)
                        power_tmp = this_trial.power_levels(i_loc);
                        for i=1:size(cell_and_pos,1)
                            if  i_cell == cell_and_pos(i,1) % Update one cell in this big for-loop
                                if strcmp(prior_info.GP_params.type,'xy_square')
                                    i_xy= cell_and_pos(i,2);i_z= cell_and_pos(i,3);
                                    if cell_and_status(i,2)==0
                                        this_xy= new_shape_samples{s}(i_cell).xy(i_xy);
                                    else
                                        this_xy= this_sample.xy(i_xy);
                                    end
                                    if cell_and_status(i,3)==0
                                        this_z= new_shape_samples{s}(i_cell).z(i_z);
                                    else
                                        this_z= this_sample.z(i_z);
                                    end
                                    
                                    stim=stim+...
                                        power_tmp*this_sample.gain*this_xy*this_z;
                                    
                                else
                                    % alternative types outdated
                                end
                            end
                        end
                    end
                end
                
                if isfield(this_sample,'delay_mean')
                    delay_params=struct;
                    delay_params.delay_mean=this_sample.delay_mean;
                    delay_params.delay_var=this_sample.delay_var;
                else
                    delay_params=struct;
                    delay_params.delay_mean=0;
                    delay_params.delay_var=1e-3;
                end
                
                [intensity_tmp] = calculate_intensities(stim,delay_params,prior_info.induced_intensity,timepoints);
                intensity_records.spike(s,:)= intensity_tmp.spike;
                intensity_records.event(s,:)= intensity_tmp.event;
                stim_records(s)= stim;
            end
            
            new_trials(i_trial).predicted.intensity.spike=[new_trials(i_trial).predicted.intensity.spike; mean(intensity_records.spike)];
            new_trials(i_trial).predicted.intensity.event=[new_trials(i_trial).predicted.intensity.event; mean(intensity_records.event)];
            new_trials(i_trial).predicted.PR=[new_trials(i_trial).predicted.PR; PR_post];
            new_trials(i_trial).predicted.source=[new_trials(i_trial).predicted.source; i_cell];
            new_trials(i_trial).predicted.stim=[ new_trials(i_trial).predicted.stim; mean(stim_records)];
        end
    end
end
%% Visualize the predicted event times
% if params.prediction
%     prs=[0.1 0.5 0.9];
%     for i_trial = 1:length(new_trials)
%         %  trials(i_trial).fitted.spike_times
%         new_trials(i_trial).fitted_summary.event_times=quantile(new_trials(i_trial).fitted.spike_times,prs);
%         new_trials(i_trial).fitted_summary.event_times=quantile(new_trials(i_trial).fitted.spike_times,prs);
%     end
%
%     % Scatter plot:
%     event_times_sum = zeros(n_trials,3);
%     true_event_time = zeros(n_trials,1);
%     for i_trial = 1:length(new_trials)
%         event_times_sum(i_trial,:)=new_trials(i_trial).fitted_summary.event_times;
%         if isempty(new_trials(i_trial).event_times)
%             true_event_time(i_trial)=Tmax;
%         else
%             true_event_time(i_trial)=new_trials(i_trial).event_times(1);
%         end
%         [unique_loc, ia,unique_indices]=unique(all_locations, 'rows');
%         %
%         figure
%         n_unq=size(unique_loc,1);
%         sub_row = 3;
%         sub_col = ceil(n_unq/sub_row);
%         [unique_pow, ~, pow_ind]=unique([trials.power_levels]);
%         color_list=lines(length(unique_pow));
%         for i_unq = 1:n_unq
%             subplot(sub_row,sub_col,i_unq)
%             these_trials = trials(unique_indices==i_unq);
%             these_indices=find(unique_indices==i_unq); % Should give trial an ID
%             % gather the power information:
%
%
%             for i_trial=1:length(these_trials) %'MarkerEdgeColor',color_list(pow_ind(idx),:),...
%                 idx=these_indices(i_trial);
%                 if ~no_spike(idx) && true_event_time(idx) > event_times_sum(idx,1) && true_event_time(idx) < event_times_sum(idx,3)
%                     scatter(true_event_time(idx)/time_factor,1.0,10,color_list(pow_ind(idx),:),'.')%event_times_sum(idx,2)/time_factor
%                 elseif ~no_spike(idx) && (true_event_time(idx) < event_times_sum(idx,1) || true_event_time(idx) > event_times_sum(idx,3))
%                     scatter(true_event_time(idx)/time_factor,2.0,10,color_list(pow_ind(idx),:),'x')%event_times_sum(idx,2)/time_factor
%                 elseif no_spike(idx) && ~(prior_info.induced_intensity.spike_time_max < event_times_sum(idx,1) || 30 > event_times_sum(idx,3))
%                     scatter(true_event_time(idx)/time_factor,2.0,10,color_list(pow_ind(idx),:),'x')%event_times_sum(idx,2)/time_factor
%                 else
%                     scatter(true_event_time(idx)/time_factor,1.0,10,color_list(pow_ind(idx),:),'.')%event_times_sum(idx,2)/time_factor
%                 end
%                 hold on;
%                 line([true_event_time(idx) true_event_time(idx)]/time_factor,event_times_sum(idx,[1 3])/time_factor,'LineStyle',':','LineWidth',1.0,'Color',color_list(pow_ind(idx),:))
%                 hold on;
%             end
%             xlim([0 Tmax]/time_factor); ylim([0 Tmax]/time_factor);
%             %     xlim([min(true_event_time(these_indices)) max(true_event_time(these_indices))+1]);
%             %    ylim([min(event_times_sum(these_indices,1)) max(event_times_sum(these_indices, 3))+1]);
%             line([0 Tmax]/time_factor, [0 Tmax]/time_factor)
%             xlabel('Observed (ms)');ylabel('Predicted (ms)')
%             title([num2str(round(unique_loc(i_unq,1)-neurons(1).location(1),1) ) ' ' num2str(round(unique_loc(i_unq,2)-neurons(1).location(2),1)) ' ' num2str(round(unique_loc(i_unq,3)-neurons(1).location(3),1))]);
%             %     for i_pow = 1:length(unique_pow)
%             %         txt_string = ['Power ' num2str(unique_pow(i_pow))];
%             %         text(Tmax*0.8/time_factor,10*i_pow/time_factor,txt_string,'Color', color_list(i_pow,:))
%             %     end
%         end
%
%     end
% end