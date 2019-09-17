function [trials] = generate_psc_data_batch(neurons,trials,stimuli_size,prior_info,simulation_params)
%% 
background_rate=simulation_params.background_rate;
spike_time_max=prior_info.induced_intensity.spike_time_max;
event_time_max=prior_info.induced_intensity.event_time_max;
delay_indicator=simulation_params.batch.delay_indicator;
%%
mu_bg =1/background_rate;
number_of_trials = length(trials);
number_of_cells =length(neurons);
temp_struct=struct;temp_struct.event_times=[];temp_struct.spike_times=[];temp_struct.assignments=[];
for i_trial = 1:number_of_trials
 trials(i_trial).truth=temp_struct;
    for i_cell = 1:number_of_cells
        k=stimuli_size(i_trial,i_cell);
        params_sim=neurons(i_cell).truth;
%         params_sim.linkfuncs=experiment_setup.prior_info.induced_intensity.linkfunc;
        stim_threshold=prior_info.induced_intensity.minimum_stim_threshold/neurons(i_cell).truth.gain;
        if k > stim_threshold
            stim=k*neurons(i_cell).truth.gain;
            [spikes,events] = spike_curves_sim(stim,params_sim,prior_info.induced_intensity);
            
            if ~isempty(spikes)
                for i = 1:length(spikes)
                    if rand(1) < params_sim.PR
                        % censoring at time max:
                        if events(i)<event_time_max & spikes(i)<spike_time_max
                            trials(i_trial).truth.event_times=[trials(i_trial).truth.event_times events(i)];
                            trials(i_trial).truth.spike_times=[trials(i_trial).truth.spike_times spikes(i)];
                            trials(i_trial).truth.assignments=[trials(i_trial).truth.assignments neurons(i_cell).cell_ID];

                        end
                    end
                end
            end
        end
    end
    
    % add background event:
    R = exprnd(mu_bg);
    while R < event_time_max
        trials(i_trial).truth.spike_times=[trials(i_trial).truth.spike_times max(1,round(R))];
        trials(i_trial).truth.event_times=[trials(i_trial).truth.event_times max(1,round(R))];
        
        trials(i_trial).truth.assignments=[trials(i_trial).truth.assignments 0];
        R = R+exprnd(mu_bg);
    end
    if delay_indicator
        trials(i_trial).event_times = trials(i_trial).truth.event_times;
    else
         trials(i_trial).event_times = trials(i_trial).truth.spike_times;
    end
    if mod(i_trial,50)==1
    fprintf('Trials generated: %d;\n',i_trial)
    end
end
