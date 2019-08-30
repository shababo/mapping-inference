function [trials ] = get_fitted_values(neurons, trials, prior_info, params)
%% Outline:
% - Draw samples of all relevant parameters from their posterior
% distributions
% - Draw spike time and event time given each sample
% - Summarize the posterior samples (parameters and event times)

% params.plot_only
% params.MC_params.sample_size = 50;
n_cell = length(neurons);
Tmax=prior_info.induced_intensity.event_time_max;
time_factor = prior_info.induced_intensity.time_factor;
n_trials = length(trials);
timepoints = 1:Tmax;
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

%% Calculate the fitted intensities (as posterior averages)
if isfield(posterior_samples,'background')
    bg_tmp =cellfun(@(x) x.background,posterior_samples);
    background_post=mean(bg_tmp);
else
    background_post = 0;
end
    

for i_trial = 1:n_trials
    trials(i_trial).fitted=struct;
    trials(i_trial).fitted.intensity=struct;
    trials(i_trial).fitted.intensity.spike= background_post*ones(1,length(timepoints));
    trials(i_trial).fitted.intensity.event=background_post*ones(1,length(timepoints));
    trials(i_trial).fitted.PR=1;
    trials(i_trial).fitted.source=0;
    trials(i_trial).fitted.timepoints=timepoints;
    trials(i_trial).fitted.time_factor=time_factor;
    trials(i_trial).fitted.stim=0;
    
end

for i_cell = 1:n_cell
    
    if isfield(posterior_samples,'PR')
        PR_tmp =cellfun(@(x) x.PR,posterior_samples);
        PR_post=mean(PR_tmp);
    else
        PR_post = 1.0;
    end
    
    for i_trial =1:n_trials
        this_trial=trials(i_trial);
        %         rel_pos=neurons(i_cell).params(end).shapes.locations - ...
        %             ones(size(neurons(i_cell).params(end).shapes.locations,1),1)*(this_trial.locations-neurons(i_cell).location);
        %         i_shape=find(sum( (rel_pos.^2)')==0);
        intensity_records=struct;
        intensity_records.spike= zeros(S,length(timepoints));
        intensity_records.event= zeros(S,length(timepoints));
        stim_records=zeros(S,1);
        for s=1:S
            this_sample=posterior_samples{s}(i_cell);
            for i_loc = 1:size(trials(i_trial).locations,1)
                cell_and_pos=trials(i_trial).cell_and_pos{i_loc};
                stim=0;
                %                 this_trial.power_levels*this_sample.gain*this_sample.shapes(i_shape);
                
                if ~isempty(cell_and_pos)
                    power_tmp = this_trial.power_levels(i_loc);
                    if ~(isfield(neurons(1).params,'shapes') | isfield(neurons(1).params,'xy'))
                        for i=1:length(cell_and_pos)
                            if  i_cell == cell_and_pos(i)
                                stim=stim+...
                                    power_tmp*this_sample.gain;
                            end
                        end
                    else
                        for i=1:size(cell_and_pos,1)
                            if  i_cell == cell_and_pos(i,1) % Update one cell in this big for-loop
                                if strcmp(prior_info.GP_params.type,'xy_square')
                                    i_xy= cell_and_pos(i,2);i_z= cell_and_pos(i,3);
                                    stim=stim+...
                                        power_tmp*this_sample.gain*this_sample.xy(i_xy)*this_sample.z(i_z);
                                else
                                    i_pos= cell_and_pos(i,2);
                                    stim=stim+...
                                        power_tmp*this_sample.gain*this_sample.shapes(i_pos);
                                end
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
        
        trials(i_trial).fitted.intensity.spike=[trials(i_trial).fitted.intensity.spike; mean(intensity_records.spike)];
        trials(i_trial).fitted.intensity.event=[trials(i_trial).fitted.intensity.event; mean(intensity_records.event)];
        trials(i_trial).fitted.PR=[trials(i_trial).fitted.PR; PR_post];
        trials(i_trial).fitted.source=[trials(i_trial).fitted.source; i_cell];
        trials(i_trial).fitted.stim=[trials(i_trial).fitted.stim; mean(stim_records)];
        
        
    end
end
