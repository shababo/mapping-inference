function [trials]=get_posterior_event_probability(neurons, trials, prior_info)
%% Outline:
% - Draw samples of all relevant parameters from their posterior
% distributions
% - Draw additional parameters given the other parameters (e.g., shape
% value at a new location
% -
n_cell = length(neurons);
Tmax=prior_info.induced_intensity.event_time_max;
time_factor = 20;
S= 100;
%% Extract posterior information from neurons
params.prediction=false;
clear('posterior_params')
for i_cell = 1:n_cell
    posterior_params(i_cell)=neurons(i_cell).params(end);
end
%% Draw samples from posterior distributions
posterior_samples = cell([S 1]);
for s =1:S
    [posterior_samples{s},~] = draw_samples_from_var_dist(posterior_params);
end
%%
spike_curves=prior_info.induced_intensity;
current_lb=min(spike_curves.current);
current_gap=spike_curves.current(2)-spike_curves.current(1);
current_max_grid = length(spike_curves.current);
spike_curves_var=spike_curves.sd.^2;
spike_curves_mean=spike_curves.mean;
pre_density=struct;
grid.bound=4;grid.gap=0.1;
pre_density.grid=grid;
pre_density.normal_grid = -grid.bound:grid.gap:grid.bound;
for i = 1:length(pre_density.normal_grid)
    pre_density.cdf_grid(i) =normcdf(pre_density.normal_grid(i),0,1);
    pre_density.pdf_grid(i)=normpdf(pre_density.normal_grid(i),0,1);
end
max_grid =length(pre_density.pdf_grid);

%% Draw spike times on the fitted trials
n_trials = length(trials);
for i_cell = 1:n_cell
    for i_trial =1:n_trials
        this_trial=trials(i_trial);
        if isempty(trials(i_trial).event_times)
            trials(i_trial).post_density=[];
        else
            event_times=trials(i_trial).event_times;
            density_events = zeros(S,length(trials(i_trial).event_times));
            for s=1:S
                this_sample=posterior_samples{s}(i_cell);
                
                for i_loc = 1:size(trials(i_trial).locations)
                    stim=0;
                    cell_and_pos=trials(i_trial).cell_and_pos{i_loc};
                    if ~isempty(cell_and_pos)
                        power_tmp = this_trial.power_levels(i_loc);
                        
                        for i=1:size(cell_and_pos,1)
                            if  i_cell == cell_and_pos(i,1) % Update one cell in this big for-loop
                                
                                if strcmp(prior_info.prior_parameters.GP_params.type,'xy_square')
                                    i_xy= cell_and_pos(i,2);i_z= cell_and_pos(i,3);
                                    stim=stim+...
                                        power_tmp*this_sample.gain*this_sample.xy(i_xy)*this_sample.z(i_z);
                                    
                                else
                                    i_pos= cell_and_pos(i,2);
                                    stim=stim+...
                                        power_tmp*this_sample.gain*this_samples.shape(i_pos);
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                if isfield(this_sample,'delay_mu')
                    delay_params=struct;
                    delay_params.delay_mean=this_sample.delay_mu;
                    delay_params.delay_var=this_sample.delay_sigma^2;
                else
                    delay_params=struct;
                    delay_params.delay_mean=0;
                    delay_params.delay_var=1e-3;
                end
                
                
                
                delay_mu_temp=delay_params.delay_mean;
                delay_var_temp= delay_params.delay_var;
                stim_index= min(current_max_grid,...
                    max(1,round((stim-current_lb)/current_gap)));
                spike_times_cond_shape=spike_curves_mean(stim_index);
                expectation=delay_mu_temp+spike_times_cond_shape;
                standard_dev=sqrt(delay_var_temp+  mean(spike_curves_var(stim_index)));
                
                
                pdf_index = max(1,min(max_grid,round( ((event_times-expectation)/standard_dev +grid.bound)/grid.gap)));
                density_events(s,:)=this_sample.PR*pre_density.pdf_grid(pdf_index)/standard_dev;
                
            end
            if ~isfield(trials(i_trial),'post_density')
                trials(i_trial).post_density=zeros(n_cell+1,length(trials(i_trial).event_times));
               
            end
            trials(i_trial).post_density(i_cell+1,:)=mean(density_events);
            
        end
         trials(i_trial).post_density(1,:)=prior_info.background_rate;
    end
end
%%
for i_trial = 1:n_trials
    total_prob =sum(trials(i_trial).post_density);
    trials(i_trial).soft_assignments =  trials(i_trial).post_density;
    if ~isempty(trials(i_trial).post_density)
        for i_source = 1:size(trials(i_trial).soft_assignments,1)
            trials(i_trial).soft_assignments(i_source,:) =  trials(i_trial).post_density(i_source,:)./total_prob;
        end
    end
end



