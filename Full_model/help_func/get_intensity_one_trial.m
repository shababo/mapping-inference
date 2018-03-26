function [trial_intensity] = get_intensity_one_trial(trial, neighbourhood,experiment_setup,varargin)

% trial, neighbourhood,experiment_setup,thresholding
% neighbourhood = neighbourhoods(2,9);
% experiment_queries(2,2)
if ~isempty(varargin)
   thresholding = varargin{1};
else
    thresholding = false;
end
spike_curves=experiment_setup.prior_info.induced_intensity;
%%

good_spots = find(~isnan(trial.cell_IDs));
stim_size = zeros(length(neighbourhood.neurons),1);
for i_spot = good_spots
   i_cell_this_hood= find([neighbourhood.neurons(:).cell_ID]==trial.cell_IDs(i_spot));
   stim_size=stim_size+trial.power_levels(i_spot)*...
       neighbourhood.neurons(i_cell_this_hood).stim_locations.(trial.group_ID).effect(:,trial.location_IDs(i_spot)); 
end

%% Calculate the actual intensity received by neurons 
batch_ID=neighbourhood.batch_ID;
neurons=neighbourhood.neurons;
properties={'PR_params','gain_params','delay_mu_params','delay_sigma_params'};summary_stat={'mean'};
temp_output=grab_values_from_neurons(batch_ID,neurons,properties,summary_stat);
gamma_mean=temp_output.PR_params.mean;
gain_mean=temp_output.gain_params.mean;
delay_mu_mean=temp_output.delay_mu_params.mean;
delay_sigma_mean=temp_output.delay_sigma_params.mean;

%%
stim_received= gain_mean.*stim_size;

minimum_stim_threshold=experiment_setup.prior_info.induced_intensity.minimum_stim_threshold/3;
if thresholding
stimulated_cells = find(stim_received>minimum_stim_threshold);
else
    stimulated_cells = find(stim_received>0);
end
stim_effective = stim_received(stimulated_cells);
gamma_effective=gamma_mean(stimulated_cells);
delay_mean_effective=delay_mu_mean(stimulated_cells);
delay_sigma_effective=delay_sigma_mean(stimulated_cells);

%%
% Use evenly-spaced grid for spike_curves.current for easy mapping:
current_lb=min(spike_curves.current);
current_gap=spike_curves.current(2)-spike_curves.current(1);
current_max_grid = length(spike_curves.current);
time_grid = 1:(2*spike_curves.time_max);
trial_intensity=struct;
trial_intensity.stim_neurons = struct;
estimated_intensity=zeros(1,length(1:length(time_grid)));
stim_index=zeros(length(stimulated_cells),1);
for i_cell = 1:length(stimulated_cells)
    trial_intensity.stim_neurons(i_cell)=struct;
end
for i_stim = 1:length(stimulated_cells)
    stim_index(i_stim)= min(current_max_grid,...
        max(1,round((stim_effective(i_stim)-current_lb)/current_gap)));
    expectation=delay_mean_effective(i_stim)+spike_curves.mean(stim_index(i_stim));
    standard_dev=sqrt(delay_sigma_effective(i_stim)^2+...
        spike_curves.sd(stim_index(i_stim))^2);
    
    i_cell_this_hood=stimulated_cells(i_stim);
    trial_intensity.stim_neurons(i_stim).cell_ID=neighbourhood.neurons(i_cell_this_hood).cell_ID;
    trial_intensity.stim_neurons(i_stim).PR=gamma_mean(i_cell_this_hood);
    trial_intensity.stim_neurons(i_stim).gain=gain_mean(i_cell_this_hood);
    trial_intensity.stim_neurons(i_stim).intensity=normpdf(time_grid,expectation,standard_dev);
    trial_intensity.stim_neurons(i_stim).stimulation=stim_effective(i_stim)/gain_mean(i_cell_this_hood); 
   estimated_intensity=estimated_intensity+ trial_intensity.stim_neurons(i_stim).intensity*trial_intensity.stim_neurons(i_stim).PR;
end
estimated_intensity=estimated_intensity+experiment_setup.patched_neuron.background_rate;
%% Total intensity

% Time-rescaled version (accumulated intensity)
scaled_times=[];
if ~isempty(trial.event_times)
    events= [0 trial.event_times];
    
    for i_event = 2:length(events)
        if trial.event_times(i_event-1) < spike_curves.time_max
            scaled_times(1)=sum(estimated_intensity( (events(i_event-1)+1):events(i_event)));
        end
    end
end
%%
trial_intensity.scaled_times=scaled_times;
trial_intensity.event_times=trial.event_times;

trial_intensity.estimated_intensity=estimated_intensity;
trial_intensity.background_rate= experiment_setup.patched_neuron.background_rate;
trial_intensity.trial_ID=trial.trial_ID;
