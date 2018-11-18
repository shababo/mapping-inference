function [trials ]=get_fitted_values(neurons, trials, new_trials,prior_info,params)
%% Outline:
% - Draw samples of all relevant parameters from their posterior
% distributions
% - Draw additional parameters given the other parameters (e.g., shape
% value at a new location
% - Draw spike time and event time given each sample
% - Summarize the posterior samples (parameters and event times)
% - Visualize the posterior samples v.s. true values
% - Visualize the predicted/fitted spike times v.s. true times

n_cell = length(neurons);
Tmax=200;%prior_info.induced_intensity.event_time_max;
time_factor = 20;
%% Extract posterior information from neurons
if isempty(new_trials)
    params.prediction=false;
else
    params.prediction=true;
end
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
%% Draw additional parameters:
if params.prediction
%     [new_shape_params]=get_new_shape_conditional(neurons,new_trials,prior_info);
            [new_shape_params]=get_new_shape_conditional(current_params,new_locations,prior_info);

            
    if size(new_shape_params,1)==0
        new_shape_samples = cell([S 1]);
        for s=1:S
            [new_shape_samples{s}]=draw_samples_from_shape_conditional(new_shape_params,posterior_samples{s});
        end
    else
        new_shape_samples=cell([S 1]);
        emp_struct=struct;
        for i_cell = 1:n_cell
           emp_struct(i_cell).shapes=[]; 
        end
        for s=1:S
            new_shape_samples{s}=emp_struct;%draw_samples_from_shape_conditional(new_shape_params,posterior_samples{s});
        end
    end
end
%% Draw spike times on the fitted trials
n_trials = length(trials);
for i_cell = 1:n_cell
    for i_trial =1:n_trials
        this_trial=trials(i_trial);
        rel_pos=neurons(i_cell).params(end).shapes.locations - ...
            ones(size(neurons(i_cell).params(end).shapes.locations,1),1)*(this_trial.locations-neurons(i_cell).location);
        i_shape=find(sum( (rel_pos.^2)')==0);
        spike_records = zeros(S,1);
        event_records=zeros(S,1);
        for s=1:S
            this_sample=posterior_samples{s}(i_cell);
            stim=this_trial.power_levels*this_sample.gain*this_sample.shapes(i_shape);
            if isfield(this_sample,'delay_mu')
                delay_params=struct;
                delay_params.delay_mean=this_sample.delay_mu;
                delay_params.delay_var=this_sample.delay_sigma^2;
                
            else
                delay_params=struct;
                delay_params.delay_mean=0;
                delay_params.delay_var=1e-3;
            end
            [spikes,events] = spike_curves_sim(stim,delay_params,prior_info.induced_intensity);
            if isempty(spikes)
                spike_records(s)=NaN;
                event_records(s)=NaN;
            else
                spike_records(s)=spikes;
                event_records(s)=events;
            end
        end
        if ~isfield(trials(i_trial),'fitted')
            trials(i_trial).fitted=struct;
            trials(i_trial).fitted.spike_times=spike_records;
            trials(i_trial).fitted.event_times=event_records;
            trials(i_trial).fitted.assignments=i_cell*ones(length(event_records),1);
        elseif isempty(trials(i_trial).fitted)
            trials(i_trial).fitted=struct;
            trials(i_trial).fitted.spike_times=spike_records;
            trials(i_trial).fitted.event_times=event_records;
            trials(i_trial).fitted.assignments=i_cell*ones(length(event_records),1);
        else
            trials(i_trial).fitted.spike_times=[trials(i_trial).fitted.spike_times; spike_records];
            trials(i_trial).fitted.event_times=[trials(i_trial).fitted.event_times; event_records];
            trials(i_trial).fitted.assignments=[trials(i_trial).fitted.assignments; i_cell*ones(length(event_records),1)];
        end
    end
end
%% Draw spike & event time on new trials 
if params.prediction
for i_cell = 1:n_cell
    merged_locations=[neurons(i_cell).params(end).shapes.locations; new_shape_params(i_cell).locations];
    for i_trial =1:length(new_trials)
        this_trial=new_trials(i_trial);
        rel_pos=merged_locations- ...
            ones(size(merged_locations,1),1)*(this_trial.locations-neurons(i_cell).location);
        i_shape=find(sum( (rel_pos.^2)')==0);
        spike_records = zeros(S,1);
        event_records=zeros(S,1);
        for s=1:S
            this_sample=posterior_samples{s}(i_cell);
            merged_shapes= [this_sample.shapes; new_shape_samples{s}(i_cell).shapes];
            stim=this_trial.power_levels*this_sample.gain*merged_shapes(i_shape);
            if isfield(this_sample,'delay_mu')
                delay_params=struct;
                delay_params.delay_mean=this_sample.delay_mu;
                delay_params.delay_var=this_sample.delay_sigma^2;
                
            else
                delay_params=struct;
                delay_params.delay_mean=0;
                delay_params.delay_var=1e-3;
                
            end
            [spikes,events] = spike_curves_sim(stim,delay_params,prior_info.induced_intensity);
            if isempty(spikes)
                disp('found empty')
                spike_records(s)=Inf;
                event_records(s)=Inf;
            else
                spike_records(s)=spikes;
                event_records(s)=events;
            end
        end
        if ~isfield(new_trials(i_trial),'fitted')
            new_trials(i_trial).fitted=struct;
            new_trials(i_trial).fitted.spike_times=spike_records;
            new_trials(i_trial).fitted.event_times=event_records;
            new_trials(i_trial).fitted.assignments=i_cell*ones(length(event_records),1);
        elseif isempty(new_trials(i_trial).fitted)
            new_trials(i_trial).fitted=struct;
            new_trials(i_trial).fitted.spike_times=spike_records;
            new_trials(i_trial).fitted.event_times=event_records;
            new_trials(i_trial).fitted.assignments=i_cell*ones(length(event_records),1);
        else
            new_trials(i_trial).fitted.spike_times=[new_trials(i_trial).fitted.spike_times; spike_records];
            new_trials(i_trial).fitted.event_times=[new_trials(i_trial).fitted.event_times; event_records];
            new_trials(i_trial).fitted.assignments=[new_trials(i_trial).fitted.assignments; i_cell*ones(length(event_records),1)];
        end
    end
end
end
%% Visualization:
%     - Fits
%         Draw scatter plot that visualize the shape estimates v.s. the true values
%         Draw quantile plots for estimates v.s. the true values
%         Draw predicted spike time (mean + quantiles) v.s. the true spike times
%         Draw histograms of fitted spike distributions v.s. the true spike times
%     - Prediction
%         Draw predicted shape values (distributions) at new locations (compared with true values if available)
%           - need to figure out the conditional distribution for this prediction
%         Draw predicted spike time distributions

% True values of shapes:
% trials(1).truth

%     params_delay.delay_mean=0; params_delay.delay_var=1; % not used
%
%     % Now find the true spike time:
%     % rel_loc =  pred_trials.locations - neurons(1).truth.location;
%     % this_size = griddata(simulation_params.mesh_grid(:,1),simulation_params.mesh_grid(:,2),simulation_params.mesh_grid(:,3),...
%     %     neurons(i_cell).truth.shape,rel_loc(1),rel_loc(2),rel_loc(3));
%     % stim=pred_trials.power_levels*this_size*neurons(1).truth.optical_gain;
%     %  [spikes,events] = spike_curves_sim(stim,params_sim,prior_info.induced_intensity);
%     true_spike=pred_trials.event_times;
%
%     % Visualize the fits
%     spike_records(spike_records>500)=500;
%     figure(1)
%     histogram(spike_records,40)
%     hold on;
% plot([true_spike true_spike], [0 100],'Color','k','LineWidth',3)
% xlabel('spike time (histogram: posterior samples; bar: true spike time)');
% xlim([0 500])
% ylabel('Frequency')
% title(['Exp ' num2str(i_exp) '; Trial ' num2str(i_trial) '; Cell ' num2str(1)])
%% Visualize predicted and true spike times 
% median spike times and quantiles v.s. true spike times
prs=[0.05 0.5 0.95];
for i_trial = 1:n_trials
    %  trials(i_trial).fitted.spike_times
    trials(i_trial).fitted_summary.spike_times=quantile(trials(i_trial).fitted.spike_times,prs);
    trials(i_trial).fitted_summary.event_times=quantile(trials(i_trial).fitted.event_times,prs);
end


%% Scatter plot:

event_times_sum = zeros(n_trials,3);
true_event_time = zeros(n_trials,1);
no_spike = zeros(n_trials,1);
for i_trial = 1:n_trials
    event_times_sum(i_trial,:)=trials(i_trial).fitted_summary.event_times;
    if isempty(trials(i_trial).event_times)
        true_event_time(i_trial)=event_times_sum(i_trial,2);
        no_spike(i_trial) = 1;
    else
        true_event_time(i_trial)=trials(i_trial).event_times(1);
    end
end

figure
% hold on;
% subplot(2,2,4)
hit_count = 0;
miss_count = 0;
fa_count = 0;
correj_count = 0;
for i_trial = 1:n_trials
%     if 
    % observed spike within CI
    if ~no_spike(i_trial) && true_event_time(i_trial) < prior_info.induced_intensity.spike_time_max
        if true_event_time(i_trial) > event_times_sum(i_trial,1) && true_event_time(i_trial) < event_times_sum(i_trial,3)
            subplot(2,4,1)
            scatter(true_event_time(i_trial)/time_factor,.25 + rand*.67,15,'MarkerEdgeColor','b','MarkerFaceColor','b')
            hold on
            line([true_event_time(i_trial) true_event_time(i_trial)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','b')
            hit_count = hit_count+1;
            subplot(2,4,5)
            scatter3(trials(i_trial).locations(1)-neurons(1).location(1),trials(i_trial).locations(2)-neurons(1).location(2),...
                -trials(i_trial).locations(3)+neurons(1).location(3),15,'MarkerEdgeColor','b','MarkerFaceColor','b')
        % observed spike outside of CI
        elseif (true_event_time(i_trial) < event_times_sum(i_trial,1) || true_event_time(i_trial) > event_times_sum(i_trial,3))
            subplot(2,4,2)
            scatter(true_event_time(i_trial)/time_factor,1.25 + rand*.67,15,'MarkerEdgeColor','r','MarkerFaceColor','r')
            hold on
            line([true_event_time(i_trial) true_event_time(i_trial)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','r')
            miss_count = miss_count+1;
            subplot(2,4,6)
            scatter3(trials(i_trial).locations(1)-neurons(1).location(1),trials(i_trial).locations(2)-neurons(1).location(2),...
                -trials(i_trial).locations(3)+neurons(1).location(3),15,'MarkerEdgeColor','r','MarkerFaceColor','r')
        end
%     elseif ~no_spike(i_trial)
%         if true_event_time(i_trial) > event_times_sum(i_trial,1) && true_event_time(i_trial) < event_times_sum(i_trial,3)
%             subplot(2,2,3)
%             scatter(true_event_time(i_trial)/time_factor,.25 + rand*.67,15,'MarkerEdgeColor','c','MarkerFaceColor','c')
%             hold on
%             line([true_event_time(i_trial) true_event_time(i_trial)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','c')
%         % observed spike outside of CI
%         else%if true_event_time(i_trial) < event_times_sum(i_trial,1) || true_event_time(i_trial) > event_times_sum(i_trial,3)
%             subplot(2,2,4)
%             scatter(true_event_time(i_trial)/time_factor,1.25 + rand*.67,15,'MarkerEdgeColor','m','MarkerFaceColor','m')
%             hold on
%             line([true_event_time(i_trial) true_event_time(i_trial)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','m')
%         end
    elseif no_spike(i_trial) 
        % no observed spike but CI is within spike time bounds 
        if prior_info.induced_intensity.spike_time_max > event_times_sum(i_trial,2)
            subplot(2,4,4)
            scatter(event_times_sum(i_trial,2)/time_factor,1.25 + rand*.67,15,'MarkerEdgeColor','m','MarkerFaceColor','m')
            hold on
            line([event_times_sum(i_trial,2) event_times_sum(i_trial,2)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','m')
            fa_count = fa_count + 1;
            subplot(2,4,8)
            scatter3(trials(i_trial).locations(1)-neurons(1).location(1),trials(i_trial).locations(2)-neurons(1).location(2),...
                -trials(i_trial).locations(3)+neurons(1).location(3),15,'MarkerEdgeColor','m','MarkerFaceColor','m')
        elseif prior_info.induced_intensity.spike_time_max < event_times_sum(i_trial,1)% && prior_info.induced_intensity.time_max > event_times_sum(i_trial,1)
            subplot(2,4,3)
            scatter(event_times_sum(i_trial,1)/time_factor,1.25 + rand*.67,15,'MarkerEdgeColor','c','MarkerFaceColor','c')
            hold on
            line([event_times_sum(i_trial,1) event_times_sum(i_trial,1)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','b')
            correj_count = correj_count + 1;
            subplot(2,4,7)
            scatter3(trials(i_trial).locations(1)-neurons(1).location(1),trials(i_trial).locations(2)-neurons(1).location(2),...
                -trials(i_trial).locations(3)+neurons(1).location(3),15,'MarkerEdgeColor','c','MarkerFaceColor','c')
        elseif prior_info.induced_intensity.time_max < event_times_sum(i_trial,1)
%             subplot(2,2,3)
%             scatter(event_times_sum(i_trial,1)/time_factor,1.25 + rand*.67,15,'MarkerEdgeColor','c','MarkerFaceColor','c')
%             hold on
%             line([event_times_sum(i_trial,1) event_times_sum(i_trial,1)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','c')
%         else
%             subplot(2,2,2)
%             scatter(true_event_time(i_trial)/time_factor,1.25 + rand*.67,15,'MarkerEdgeColor','c','MarkerFaceColor','c')
%             hold on
%             line([true_event_time(i_trial) true_event_time(i_trial)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','c')
        end
       
    end
hold on;
end

% for ii = 1:4
    subplot(2,4,1)
%     xlim([0 Tmax]/time_factor); ylim([0 Tmax]/time_factor);
%     xlim([min(true_event_time(these_indices)) max(true_event_time(these_indices))+1]);
%    ylim([min(event_times_sum(these_indices,1)) max(event_times_sum(these_indices, 3))+1]);
     line([0 Tmax]/time_factor, [0 Tmax]/time_factor)
     line(prior_info.induced_intensity.spike_time_max*[1 1]/time_factor, [0 10])
    xlabel('Observed (ms)');ylabel('Predicted (ms)')
    
    subplot(2,4,3)
%     xlim([0 Tmax]/time_factor); ylim([0 Tmax]/time_factor);
%     xlim([min(true_event_time(these_indices)) max(true_event_time(these_indices))+1]);
%    ylim([min(event_times_sum(these_indices,1)) max(event_times_sum(these_indices, 3))+1]);
     line([0 Tmax]/time_factor, [0 Tmax]/time_factor)
     line(prior_info.induced_intensity.spike_time_max*[1 1]/time_factor, [0 10])
    xlabel('Observed (ms)');ylabel('Predicted (ms)')
    
    subplot(2,4,2)
%      ylim([0 Tmax]/time_factor);
%     xlim([min(true_event_time(these_indices)) max(true_event_time(these_indices))+1]);
%    ylim([min(event_times_sum(these_indices,1)) max(event_times_sum(these_indices, 3))+1]);
     line([0 Tmax]/time_factor, [0 Tmax]/time_factor)
%      line([0 20],[0 20])
    xlabel('Observed (ms)');ylabel('Predicted (ms)')
    
    subplot(2,4,4)
%     xlim([0 Tmax]/time_factor); ylim([0 Tmax]/time_factor);
%     xlim([min(true_event_time(these_indices)) max(true_event_time(these_indices))+1]);
%    ylim([min(event_times_sum(these_indices,1)) max(event_times_sum(these_indices, 3))+1]);
     line([0 Tmax]/time_factor, [0 Tmax]/time_factor)
     line(prior_info.induced_intensity.spike_time_max*[1 1]/time_factor, [0 10])
    xlabel('Observed (ms)');ylabel('Predicted (ms)')
% end

hit_count
miss_count
fa_count
correj_count

% figure
% % hold on;
% % subplot(2,2,4)
% for i_trial = 1:n_trials
% %     if 
%     % observed spike within CI
%     if ~no_spike(i_trial) && true_event_time(i_trial) > event_times_sum(i_trial,1) && true_event_time(i_trial) < event_times_sum(i_trial,3)
%         scatter3(trials(i_trial).locations(1),trials(i_trial).locations(2),trials(i_trial).locations(3),15,'MarkerEdgeColor','b','MarkerFaceColor','b')
%         hold on
%     % observed spike outside of CI
%     elseif ~no_spike(i_trial)  && (true_event_time(i_trial) < event_times_sum(i_trial,1) || true_event_time(i_trial) > event_times_sum(i_trial,3))
%         scatter3(trials(i_trial).locations(1),trials(i_trial).locations(2),trials(i_trial).locations(3),15,'MarkerEdgeColor','r','MarkerFaceColor','r')
%         hold on
%     % no observed spike but CI is within spike time bounds    
%     elseif no_spike(i_trial) && ~(prior_info.induced_intensity.spike_time_max < event_times_sum(i_trial,1) || 30 > event_times_sum(i_trial,3))
%         scatter3(trials(i_trial).locations(1),trials(i_trial).locations(2),trials(i_trial).locations(3),15,'MarkerEdgeColor','m','MarkerFaceColor','m')
%         hold on
% %         line([true_event_time(i_trial) true_event_time(i_trial)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','m')
%     elseif no_spike(i_trial) && ~(200 < event_times_sum(i_trial,1))
%         scatter3(trials(i_trial).locations(1),trials(i_trial).locations(2),trials(i_trial).locations(3),15,'MarkerEdgeColor','c','MarkerFaceColor','c')
%         hold on
% %         line([true_event_time(i_trial) true_event_time(i_trial)]/time_factor,event_times_sum(i_trial,[1 3])/time_factor,'LineStyle',':','LineWidth',.75,'color','c')
%     end
% hold on;
% end
% 
% axis image
% %  xlim([0 Tmax]/time_factor);ylim([0 Tmax]/time_factor);
%  
%  xlabel('x (vert)');ylabel('y (horiz)');zlabel('z (axial/horiz)')
% %  title('Event time');
% % title(['Hit Count: ' 

return
%% Summarize the unique locations we stimulated:
%
% return
all_locations=zeros(0,3);
for i_trial = 1:length(trials)
    if ~isempty(trials(i_trial).event_times)
        all_locations=[all_locations; trials(i_trial).locations];
    end
end
[unique_loc, ia,unique_indices]=unique(all_locations, 'rows');
%
figure
n_unq=size(unique_loc,1);
sub_row = 3;
sub_col = ceil(n_unq/sub_row);
[unique_pow, ~, pow_ind]=unique([trials.power_levels]);
color_list=lines(length(unique_pow));
for i_unq = 1:n_unq
    subplot(sub_row,sub_col,i_unq)
    these_trials = trials(unique_indices==i_unq);
    these_indices=find(unique_indices==i_unq); % Should give trial an ID
    % gather the power information:
    
    
    for i_trial=1:length(these_trials) %'MarkerEdgeColor',color_list(pow_ind(idx),:),...
        idx=these_indices(i_trial);
        if ~no_spike(idx) && true_event_time(idx) > event_times_sum(idx,1) && true_event_time(idx) < event_times_sum(idx,3)
            scatter(true_event_time(idx)/time_factor,1.0,10,color_list(pow_ind(idx),:),'.')%event_times_sum(idx,2)/time_factor
        elseif ~no_spike(idx) && (true_event_time(idx) < event_times_sum(idx,1) || true_event_time(idx) > event_times_sum(idx,3))
            scatter(true_event_time(idx)/time_factor,2.0,10,color_list(pow_ind(idx),:),'x')%event_times_sum(idx,2)/time_factor
        elseif no_spike(idx) && ~(prior_info.induced_intensity.spike_time_max < event_times_sum(idx,1) || 30 > event_times_sum(idx,3))
            scatter(true_event_time(idx)/time_factor,2.0,10,color_list(pow_ind(idx),:),'x')%event_times_sum(idx,2)/time_factor
        else
            scatter(true_event_time(idx)/time_factor,1.0,10,color_list(pow_ind(idx),:),'.')%event_times_sum(idx,2)/time_factor
        end
        hold on;
        line([true_event_time(idx) true_event_time(idx)]/time_factor,event_times_sum(idx,[1 3])/time_factor,'LineStyle',':','LineWidth',1.0,'Color',color_list(pow_ind(idx),:))
        hold on;
    end
    xlim([0 Tmax]/time_factor); ylim([0 Tmax]/time_factor);
%     xlim([min(true_event_time(these_indices)) max(true_event_time(these_indices))+1]);
%    ylim([min(event_times_sum(these_indices,1)) max(event_times_sum(these_indices, 3))+1]);
     line([0 Tmax]/time_factor, [0 Tmax]/time_factor)
    xlabel('Observed (ms)');ylabel('Predicted (ms)')
    title([num2str(round(unique_loc(i_unq,1)-neurons(1).location(1),1) ) ' ' num2str(round(unique_loc(i_unq,2)-neurons(1).location(2),1)) ' ' num2str(round(unique_loc(i_unq,3)-neurons(1).location(3),1))]);
%     for i_pow = 1:length(unique_pow)
%         txt_string = ['Power ' num2str(unique_pow(i_pow))];
%         text(Tmax*0.8/time_factor,10*i_pow/time_factor,txt_string,'Color', color_list(i_pow,:))
%     end
end

 %% Visualize the predicted event times 
 if params.prediction
 prs=[0.1 0.5 0.9];
for i_trial = 1:length(new_trials)
    %  trials(i_trial).fitted.spike_times
    new_trials(i_trial).fitted_summary.event_times=quantile(new_trials(i_trial).fitted.spike_times,prs);
    new_trials(i_trial).fitted_summary.event_times=quantile(new_trials(i_trial).fitted.spike_times,prs);
end

% Scatter plot:
event_times_sum = zeros(n_trials,3);
true_event_time = zeros(n_trials,1);
for i_trial = 1:length(new_trials)
    event_times_sum(i_trial,:)=new_trials(i_trial).fitted_summary.event_times;
    if isempty(new_trials(i_trial).event_times)
        true_event_time(i_trial)=Tmax;
    else
        true_event_time(i_trial)=new_trials(i_trial).event_times(1);
    end
end
figure
scatter(true_event_time,event_times_sum(:,2),25,'MarkerEdgeColor','b','MarkerFaceColor','b')
hold on;
for i_trial = 1:length(new_trials)
line([true_event_time(i_trial) true_event_time(i_trial)],event_times_sum(i_trial,[1 3]),'LineStyle',':','LineWidth',0.2,'Color','b')
hold on;
end
 xlim([0 Tmax]);ylim([0 Tmax]);
 line([0 Tmax], [0 Tmax])
 xlabel('Observed (ms)');ylabel('Predicted (ms)')
 title('Event time');
 
%% Summarize the unique locations in the new trials:
% 

all_locations=zeros(0,3);
for i_trial = 1:length(new_trials)
    all_locations=[all_locations; new_trials(i_trial).locations];
end
[unique_loc, ia,unique_indices]=unique(all_locations, 'rows');
figure
n_unq=size(unique_loc,1);

for i_unq = 1:n_unq
    subplot(1,n_unq,i_unq)
    these_trials = new_trials(unique_indices==i_unq);
    these_indices=find(unique_indices==i_unq); % Should give trial an ID
    % gather the power information:
    [unique_pow, ~, pow_ind]=unique([these_trials.power_levels]);
    color_list=lines(length(unique_pow));
    for i_trial=1:length(these_trials)
        idx=these_indices(i_trial);
        scatter(true_event_time(idx),event_times_sum(idx,2),25,'MarkerEdgeColor',color_list(pow_ind(i_trial),:),...
            'MarkerFaceColor',color_list(pow_ind(i_trial),:))
        hold on;
        line([true_event_time(idx) true_event_time(idx)],event_times_sum(idx,[1 3]),'LineStyle',':','LineWidth',0.2,'Color',color_list(pow_ind(i),:))
        hold on;
        
    end
    xlim([min(true_event_time(these_indices)) max(true_event_time(these_indices))+1]);ylim([min(event_times_sum(these_indices,1)) max(event_times_sum(these_indices, 3))+1]);
    line([0 Tmax], [0 Tmax])
    xlabel('Observed (ms)');ylabel('Predicted (ms)')
    title(['Event time at ' num2str(round(unique_loc(i_unq,1),1)) ' '...
        num2str(round(unique_loc(i_unq,2),1)) ' ' num2str(round(unique_loc(i_unq,3),1))]);
    for i_pow = 1:length(unique_pow)
        txt_string = ['Power ' num2str(unique_pow(i_pow))];
        text(Tmax*0.8,30*i_pow,txt_string,'Color', color_list(i_pow,:))
    end
end
 end
%% Visualize the posteriors of shapes v.s. true shapes, and gain*shapes versus the truth
% true_shapes=cell([n_cell 1]);
% for i_cell = 1:n_cell
%     rel_loc=neurons(i_cell).params(end).shapes.locations;
%     true_shapes{i_cell} = griddata(simulation_params.mesh_grid(:,1),simulation_params.mesh_grid(:,2),simulation_params.mesh_grid(:,3),...
%     neurons(i_cell).truth.shape,rel_loc(:,1),rel_loc(:,2),rel_loc(:,3),'linear');
% end   
% 
% % Reformat the samples for shapes, and gains
% shape_samples=cell([n_cell 1]);
% shape_gain_samples=cell([n_cell 1]);
% for i_cell = 1:n_cell
%     shape_samples{i_cell}= zeros(S,length(true_shapes{i_cell}));
%     shape_gain_samples{i_cell}= zeros(S,length(true_shapes{i_cell}));
%     for s = 1:S
%         shape_samples{i_cell}(s,:)=posterior_samples{s}(i_cell).shapes;
%         shape_gain_samples{i_cell}(s,:)=posterior_samples{s}(i_cell).shapes*posterior_samples{s}(i_cell).gain;
%     end
%     
% end
% 
% % Obtain the posterior quantiles 
% shape_summary=cell([n_cell 1]);
% shape_gain_summary=cell([n_cell 1]);
% prs=[0.1 0.5 0.9];
% for i_cell = 1:n_cell
%     shape_summary{i_cell}= zeros(length(true_shapes{i_cell}),3);
%     shape_gain_summary{i_cell}= zeros(length(true_shapes{i_cell}),3);
%     for i_shape = 1:length(true_shapes{i_cell})
%     
%     shape_summary{i_cell}(i_shape,:)=quantile(shape_samples{i_cell}(:,i_shape),prs);
%     shape_gain_summary{i_cell}(i_shape,:)=quantile(shape_gain_samples{i_cell}(:,i_shape),prs);
%     end
% end
% 
% 
% figure(3)
% for i_cell =1:n_cell
% subplot(1,n_cell,i_cell)
% scatter(true_shapes{i_cell},shape_summary{i_cell}(:,2),25,'MarkerEdgeColor','b','MarkerFaceColor','b')
% hold on;
% 
% for i_shape = 1:length(true_shapes{i_cell})
% line([true_shapes{i_cell}(i_shape) true_shapes{i_cell}(i_shape)],...
%     shape_summary{i_cell}(i_shape,[1 3]),'LineStyle',':','LineWidth',0.2,'Color',color_list(pow_ind(i),:))
% hold on;
% end
% ts_range = [min(true_shapes{i_cell}) max(true_shapes{i_cell})];
%  xlim(ts_range);ylim(ts_range);
%  line(ts_range, ts_range)
%  xlabel('True values');ylabel('Posteriors')
%  title(['Shape values for fitted trials; Neuron ' num2str(i_cell)]);
% end
% 
% 
% 
% 
% figure(4)
% for i_cell =1:n_cell
% subplot(1,n_cell,i_cell)
% scatter(true_shapes{i_cell}*neurons(i_cell).truth.gain,shape_gain_summary{i_cell}(:,2),25,'MarkerEdgeColor','b','MarkerFaceColor','b')
% hold on;
% for i_shape = 1:length(true_shapes{i_cell})
% line([true_shapes{i_cell}(i_shape)*neurons(i_cell).truth.gain true_shapes{i_cell}(i_shape)*neurons(i_cell).truth.gain],...
%     shape_gain_summary{i_cell}(i_shape,[1 3]),'LineStyle',':','LineWidth',0.2,'Color',color_list(pow_ind(i),:))
% hold on;
% end
% ts_range = [min(true_shapes{i_cell}) max(true_shapes{i_cell})]*neurons(i_cell).truth.gain;
% xlim(ts_range);ylim(ts_range);
%  line(ts_range, ts_range)
%  xlabel('True values');ylabel('Posteriors')
%  title(['Shape*gain values for fitted trials; Neuron ' num2str(i_cell)]);
% end
% %% Visualize predicted posteriors of shapes v.s. true shapes for new trials 
%   
% true_shapes=cell([n_cell 1]);
% for i_cell = 1:n_cell
%     rel_loc=new_shape_params(i_cell).locations;
%     true_shapes{i_cell} = griddata(simulation_params.mesh_grid(:,1),simulation_params.mesh_grid(:,2),simulation_params.mesh_grid(:,3),...
%     neurons(i_cell).truth.shape,rel_loc(:,1),rel_loc(:,2),rel_loc(:,3),'linear');
% end   
% 
% % Reformat the samples for shapes, and gains
% shape_samples=cell([n_cell 1]);
% shape_gain_samples=cell([n_cell 1]);
% for i_cell = 1:n_cell
%     shape_samples{i_cell}= zeros(S,length(true_shapes{i_cell}));
%     shape_gain_samples{i_cell}= zeros(S,length(true_shapes{i_cell}));
%     for s = 1:S
%         shape_samples{i_cell}(s,:)=new_shape_samples{s}(i_cell).shapes;
%         shape_gain_samples{i_cell}(s,:)=new_shape_samples{s}(i_cell).shapes*posterior_samples{s}(i_cell).gain;
%     end
%     
% end
% 
% % Obtain the posterior quantiles 
% shape_summary=cell([n_cell 1]);
% shape_gain_summary=cell([n_cell 1]);
% prs=[0.1 0.5 0.9];
% for i_cell = 1:n_cell
%     shape_summary{i_cell}= zeros(length(true_shapes{i_cell}),3);
%     shape_gain_summary{i_cell}= zeros(length(true_shapes{i_cell}),3);
%     for i_shape = 1:length(true_shapes{i_cell})
%     
%     shape_summary{i_cell}(i_shape,:)=quantile(shape_samples{i_cell}(:,i_shape),prs);
%     shape_gain_summary{i_cell}(i_shape,:)=quantile(shape_gain_samples{i_cell}(:,i_shape),prs);
%     end
% end
% 
% 
% figure(5)
% for i_cell =1:n_cell
% subplot(1,n_cell,i_cell)
% scatter(true_shapes{i_cell},shape_summary{i_cell}(:,2),25,'MarkerEdgeColor','b','MarkerFaceColor','b')
% hold on;
% for i_shape = 1:length(true_shapes{i_cell})
% line([true_shapes{i_cell}(i_shape) true_shapes{i_cell}(i_shape)],...
%     shape_summary{i_cell}(i_shape,[1 3]),'LineStyle',':','LineWidth',0.2,'Color',color_list(pow_ind(i),:))
% hold on;
% end
% ts_range = [min(true_shapes{i_cell})-0.1 max(true_shapes{i_cell})+0.1];
%  xlim(ts_range);ylim(ts_range);
%  line(ts_range, ts_range)
%  xlabel('True values');ylabel('Posteriors')
%  title(['Shape values for new trials; Neuron ' num2str(i_cell)]);
% end
% 
% 
% 
% 
% figure(6)
% for i_cell =1:n_cell
% subplot(1,n_cell,i_cell)
% scatter(true_shapes{i_cell}*neurons(i_cell).truth.gain,shape_gain_summary{i_cell}(:,2),25,'MarkerEdgeColor','b','MarkerFaceColor','b')
% hold on;
% for i_shape = 1:length(true_shapes{i_cell})
% line([true_shapes{i_cell}(i_shape)*neurons(i_cell).truth.gain true_shapes{i_cell}(i_shape)*neurons(i_cell).truth.gain],...
%     shape_gain_summary{i_cell}(i_shape,[1 3]),'LineStyle',':','LineWidth',0.2,'Color',color_list(pow_ind(i),:))
% hold on;
% end
% 
% ts_range = [min(true_shapes{i_cell})-0.1 max(true_shapes{i_cell})+0.1]*neurons(i_cell).truth.gain;
% xlim(ts_range);ylim(ts_range);
%  line(ts_range, ts_range)
%  
%  xlabel('True values');ylabel('Posteriors')
%  title(['Shape*gain values for new trials; Neuron ' num2str(i_cell)]);
% end







