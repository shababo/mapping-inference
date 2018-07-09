function [loglklh] = update_likelihood(trials,variational_samples, ...
    GP_samples, background_rate,lklh_func,spike_curves,neurons,boundary)
%% 

% tic;
% tstart=toc;

n_shape = size(GP_samples.full.samples,2);
n_cell=length(variational_samples);
n_trial=length(trials);

stim_size=zeros(n_trial,n_cell,n_shape);
for i_cell = 1:n_cell
    shifts=[variational_samples(i_cell).shift_x];
    adjusted_location(i_cell)=neurons(i_cell).location-shifts;
    relevant_trials=[];
    relevant_stim_locations=zeros(0,1);
    trial_and_stim=zeros(0,2);
    
    this_adjusted_loc=adjusted_location(i_cell);
    for i_trial = 1:length(trials)
        relevant_flag = false;
        for i_loc = 1:length(trials(i_trial).locations)
            rel_position=trials(i_trial).locations(i_loc)-this_adjusted_loc;
            this_relevant_flag = (abs(rel_position)<boundary);
            if  this_relevant_flag
                relevant_stim_locations=[relevant_stim_locations; rel_position];
                trial_and_stim=[trial_and_stim;  [i_trial i_loc]];
                relevant_flag =true;
            end
        end
        if relevant_flag
            relevant_trials = [relevant_trials i_trial];
        end
    end
    
    for i_shape = 1:n_shape
        this_shape = GP_samples.full.samples(:,i_shape);
        switch size(relevant_stim_locations,2)
            case 1
                shape_val = interp1(GP_samples.full.locations,this_shape,relevant_stim_locations);
            case 2
                shape_val = griddata(GP_samples.full.locations(:,1),GP_samples.full.locations(:,2),...
                    this_shape,relevant_stim_locations(:,1),relevant_stim_locations(:,2));
            case 3
                shape_val = griddata(GP_samples.full.locations(:,1),GP_samples.full.locations(:,2),GP_samples.full.locations(:,3),...
                    this_shape,relevant_stim_locations(:,1),relevant_stim_locations(:,2),relevant_stim_locations(:,3));
        end
         
        for i_rel = 1:length(shape_val)
            i_trial = trial_and_stim(i_rel,1);
            i_loc = trial_and_stim(i_rel,2);
            if length(trials(i_trial).power_levels)==1
                power_tmp =  trials(i_trial).power_levels;
            else
                power_tmp = trials(i_trial).power_levels(i_loc);
            end
            stim_size(i_trial,i_cell,i_shape)=stim_size(i_trial,i_cell,i_shape)+...
                power_tmp*shape_val(i_rel);
        end
    end
end
%
% t2=toc;

Tmax=spike_curves.time_max;
minimum_stim_threshold=spike_curves.minimum_stim_threshold;
% Calculate a grid for standard normal r.v. for quick CDF calculation:
grid.bound=4;grid.gap=0.1;
normal_grid = -grid.bound:grid.gap:grid.bound;
for i = 1:length(normal_grid)
    cdf_grid(i) = normcdf(normal_grid(i),0,1);
    pdf_grid(i)=normpdf(normal_grid(i),0,1);
end
loglklh=zeros(n_cell,1);
% Use evenly-spaced grid for spike_curves.current for easy mapping:
current_lb=min(spike_curves.current);
current_gap=spike_curves.current(2)-spike_curves.current(1);
current_max_grid = length(spike_curves.current);

gain_sample=reshape([variational_samples(:).gain], [n_cell 1]);
PR_sample=reshape([variational_samples(:).PR], [n_cell 1]);
delay_mu_sample=reshape([variational_samples(:).delay_mu], [n_cell 1]);
delay_sigma_sample=reshape([variational_samples(:).delay_sigma], [n_cell 1]);
%
loglklh_vec = zeros(n_trial,1);
for  i_trial = 1:n_trial
    event_times=trials(i_trial).event_times;
prob_this_trial(1,:)=background_rate*ones(1,length(event_times)+1);
        prob_this_trial(1,end)=background_rate*Tmax;
    
    for i_cell = 1:n_cell % can reduce to cell with sufficiently large stimuliaton 
        delay_mu_temp=delay_mu_sample(i_cell);
        delay_sigma_temp=delay_sigma_sample(i_cell);
        stim_index = ones(n_shape,1);
        for i_shape = 1:n_shape 
            stim_temp =stim_size(i_trial,i_cell,i_shape);
            effective_stim=stim_temp*gain_sample(i_cell);
            if effective_stim>minimum_stim_threshold
                stim_index(i_shape)= min(current_max_grid,...
                        max(1,round((effective_stim-current_lb)/current_gap)));
            end
        end
        spike_times_cond_shape=spike_curves.mean(stim_index);
        expectation=delay_mu_temp+mean(spike_times_cond_shape); % mean expectation
        standard_dev=sqrt(delay_sigma_temp^2+...
            mean(spike_curves.sd(stim_index).^2)+ var(spike_times_cond_shape));
        cdf_index = max(1,min(length(cdf_grid),round( ((Tmax-expectation)/standard_dev +grid.bound)/grid.gap)));
        pdf_index = max(1,min(length(pdf_grid),round( ((event_times-expectation)/standard_dev +grid.bound)/grid.gap)));
        prob_this_trial(i_cell,:)=...
            PR_sample(i_cell)*[pdf_grid(pdf_index) cdf_grid(cdf_index)];
    end
   
    loglklh_vec(i_trial)=  lklh_func(trials(i_trial),prob_this_trial);
    %         te=toc;time_rec(4)=time_rec(4)+te-ts;
end
loglklh=sum(loglklh_vec);
% tend=toc;
