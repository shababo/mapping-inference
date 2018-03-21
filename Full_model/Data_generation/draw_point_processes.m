function experiment_query_this_group = draw_point_processes(experiment_query_this_group, ...
        this_neighbourhood,experiment_setup)
    
time_max=length(experiment_setup.prior_info.current_template);
number_of_cells=length(this_neighbourhood.neurons);
cell_ID_list = [this_neighbourhood.neurons.cell_ID];
number_of_trials=length(experiment_query_this_group.trials);
% number_of_spots=length(experiment_query_this_group.trials(1).location_IDs);
% stimuli_size=zeros(number_of_trials,number_of_cells);
group_ID=experiment_query_this_group.trials(1).group_ID;

stimuli_size = get_stim_size(group_ID,experiment_query_this_group.trials,...
    this_neighbourhood);
% for l = 1:number_of_trials
%     for m = 1:number_of_spots
%         this_loc_ID=experiment_query_this_group.trials(l).location_IDs(m);
%         if isnan(this_loc_ID)
%         else
%             this_loc_power = experiment_query_this_group.trials(l).power_levels(m);
%             cell_ID=experiment_query_this_group.trials(l).cell_IDs(m);
%             stimuli_size(l,:) = stimuli_size(l,:)+...
%                 (this_neighbourhood.neurons(cell_ID_list == cell_ID).stim_locations.(group_ID).effect(:,this_loc_ID)*this_loc_power)';
%         end
%     end
% end

mu_bg = 1/experiment_setup.patched_neuron.background_rate;

for i_trial = 1:number_of_trials

    experiment_query_this_group.trials(i_trial).truth.event_times=[];
    experiment_query_this_group.trials(i_trial).truth.assignments=[];

    for i_cell = 1:number_of_cells
        k=stimuli_size(i_trial,i_cell);
        params_sim=this_neighbourhood.neurons(i_cell).truth;
%         params_sim.linkfuncs=experiment_setup.prior_info.induced_intensity.linkfunc;
        stim_threshold=experiment_setup.prior_info.induced_intensity.minimum_stim_threshold/experiment_setup.neurons(i_cell).truth.optical_gain;
        if k > stim_threshold

%             stim = experiment_setup.prior_info.current_template*k;
%             [~, spikes]  =lif_glm_sim(stim,params_sim);
            stim=k*experiment_setup.neurons(i_cell).truth.optical_gain;
            [spikes,events] = spike_curves_sim(stim,params_sim,experiment_setup.prior_info.induced_intensity);
            
            %         plot(V_vect)
            if ~isempty(spikes)>0
%                 switch experiment_setup.prior_info.delay.type
%                     case'gamma'
%                     shape=(experiment_setup.prior_info.delay.mean^2)/...
%                         (experiment_setup.prior_info.delay.std^2);
%                     scale = experiment_setup.prior_info.delay.mean/shape;
%                     delay_vec=round(gamrnd(shape,scale,[sum(spikes) 1]));
%                 end
%                 spikes_delay =find(spikes)+delay_vec';
                for i = 1:length(spikes)
                    if rand(1) < params_sim.PR
                        % censoring at time max:
                        if events(i)<time_max
                            experiment_query_this_group.trials(i_trial).truth.event_times=...
               [experiment_query_this_group.trials(i_trial).truth.event_times events(i)];
                            experiment_query_this_group.trials(i_trial).truth.spike_times=...
               [experiment_query_this_group.trials(i_trial).truth.event_times spikes(i)];
                            experiment_query_this_group.trials(i_trial).truth.assignments=...
[experiment_query_this_group.trials(i_trial).truth.assignments this_neighbourhood.neurons(i_cell).cell_ID];

                        end
                    end

                end
            end
        end
    end
    
    % add background event:
    R = exprnd(mu_bg);
    while R < time_max
        experiment_query_this_group.trials(i_trial).truth.event_times=[experiment_query_this_group.trials(i_trial).truth.event_times max(1,round(R))];
        experiment_query_this_group.trials(i_trial).truth.spike_times=[experiment_query_this_group.trials(i_trial).truth.event_times max(1,round(R))];
        experiment_query_this_group.trials(i_trial).truth.assignments=[experiment_query_this_group.trials(i_trial).truth.assignments 0];
        R = R+exprnd(mu_bg);
    end

end


