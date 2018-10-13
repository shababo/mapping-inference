function [parameter_path]=get_fitted_values(neurons, trials, new_trials,prior_info,params)
%% Outline:
% - Draw samples of all relevant parameters from their posterior
% distributions
% - Draw additional parameters given the other parameters (e.g., shape
% value at a new location
% - Draw spike time and event time given each sample
% - Summarize the posterior samples (parameters and event times)
% - Visualize the posterior samples v.s. true values
% - Visualize the predicted/fitted spike times v.s. true times
params.MC_params.sample_size=100;
n_cell = length(neurons);
%% Extract posterior information from neurons

clear('posterior_params')
for i_cell = 1:n_cell
    posterior_params(i_cell)=neurons(i_cell).params(end);
end

%% Draw samples from posterior distributions
S= params.MC_params.sample_size;
posterior_samples = cell([S 1]);
for s =1:S
    [posterior_samples{s},~] = draw_samples_from_var_dist(variational_params);
end
%% Draw additional parameters:
if ~isempty(new_trials)
    [new_shape_params]=get_new_shape_conditional(neurons,new_trials,prior_info);
    new_shape_samples = cell([S 1]);
    for s=1:S
        [new_shape_samples{s}]=draw_samples_from_shape_conditional(new_shape_params,posterior_samples{s});
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
            trials(i_trial).fitted.spike_times=[trials(i_trial).fitted.spike_times spike_records];
            trials(i_trial).fitted.event_times=[trials(i_trial).fitted.event_times event_records];
            trials(i_trial).fitted.assignments=[trials(i_trial).fitted.assignments i_cell*ones(length(event_records),1)];
        end
    end
end
%% Draw spike & event time on new trials 

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
                spike_records(s)=NaN;
                event_records(s)=NaN;
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
            new_trials(i_trial).fitted.spike_times=[new_trials(i_trial).fitted.spike_times spike_records];
            new_trials(i_trial).fitted.event_times=[new_trials(i_trial).fitted.event_times event_records];
            new_trials(i_trial).fitted.assignments=[new_trials(i_trial).fitted.assignments i_cell*ones(length(event_records),1)];
        end
    end
end

%%
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