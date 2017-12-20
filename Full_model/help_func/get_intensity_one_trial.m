function [trial_intensity] = get_intensity_one_trial(trial, neighbourhood,experiment_setup)

% neighbourhood = neighbourhoods(2,9);
% experiment_queries(2,2)

%%

num_spots = sum(~isnan(trial.cell_IDs));
stim_size = zeros(length(neighbourhood.neurons),1);
for i_spot = 1:num_spots
   i_cell_this_hood= find([neighbourhood.neurons(:).cell_ID]==trial.cell_IDs(i_spot));
   stim_size=stim_size+trial.power_levels(i_spot)*...
       neighbourhood.neurons(i_cell_this_hood).stim_locations.(trial.group_ID).effect(:,trial.location_IDs(i_spot)); 
end

%% Calculate the actual intensity received by neurons 
batch_ID=neighbourhood.batch_ID;
neurons=neighbourhood.neurons;
properties={'PR_params','gain_params'};summary_stat={'mean'};
temp_output=grab_values_from_neurons(batch_ID,neurons,properties,summary_stat);
gamma_mean=temp_output.PR_params.mean;
gain_mean=temp_output.gain_params.mean;
%%
stim_received= gain_mean.*stim_size;
minimum_stim_threshold=experiment_setup.prior_info.induced_intensity.minimum_stim_threshold;
stimulated_cells = find(stim_received>minimum_stim_threshold);
stim_effective = stim_received(stimulated_cells);
gamma_effective=gamma_mean(stimulated_cells);

stim_scale=experiment_setup.prior_info.induced_intensity.stim_scale;
intensity_grid=experiment_setup.prior_info.induced_intensity.intensity_grid;
stim_index=min(size(intensity_grid,1),max(1,round(stim_effective*stim_scale)));
%%
trial_intensity=struct;
trial_intensity.stim_neurons = struct;
if ~isempty(stimulated_cells)
    for i_cell = 1:length(stimulated_cells)
       trial_intensity.stim_neurons(i_cell)=struct;
    end
    % initialize the struct array before modifying them
    for i_cell = 1:length(stimulated_cells)
        i_cell_this_hood=stimulated_cells(i_cell);
       trial_intensity.stim_neurons(i_cell).cell_ID=neighbourhood.neurons(i_cell_this_hood).cell_ID;
       trial_intensity.stim_neurons(i_cell).PR=gamma_mean(i_cell_this_hood);
       trial_intensity.stim_neurons(i_cell).gain=gain_mean(i_cell_this_hood);
       trial_intensity.stim_neurons(i_cell).intensity=intensity_grid(stim_index(i_cell),:);
       trial_intensity.stim_neurons(i_cell).stimulation=stim_effective(i_cell)/gain_mean(i_cell_this_hood);
    end
end
%% Total intensity

ps_intensity =intensity_grid(stim_index,:).*(gamma_effective*ones(1,size(intensity_grid(stim_index,:),2)));
estimated_intensity= sum(ps_intensity,1)+experiment_setup.patched_neuron.background_rate;

% Time-rescaled version (accumulated intensity)
scaled_times=[];
if ~isempty(trial.event_times)
    events= [0 trial.event_times];
    for i_event = 2:length(events)
        scaled_times(1)=sum(estimated_intensity( (events(i_event-1)+1):events(i_event))); 
    end
end
%%
trial_intensity.scaled_times=scaled_times;
trial_intensity.event_times=trial.event_times;

trial_intensity.estimated_intensity=estimated_intensity;
trial_intensity.background_rate= experiment_setup.patched_neuron.background_rate;
trial_intensity.trial_ID=trial.trial_ID;
