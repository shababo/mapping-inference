function experiment_query_this_group = draw_point_processes(experiment_query_this_group, ...
        this_neighbourhood,experiment_setup)

    number_of_cells=length(experiment_setup.neurons);
    number_of_trials=length(experiment_query_this_group.trial);
    stimuli_size=zeros(number_of_trials,number_of_cells);
    this_group=experiment_query_this_group.trials(l).group_ID;
    for l = 1:number_of_trials
        for m = 1:number_of_cells
            this_loc_ID=experiment_query_this_group.trials(l).location_ID(m);
            if isnan(this_loc_ID)
            else
                this_loc_power = experiment_query_this_group.trials(l).power_levels(m);
                this_cell=experiment_query_this_group.trials(l).target_cell_ID(m);
                stimuli_size(l,:) = stimuli_size(l,:)+...
                    (this_neighbourhood.neurons(this_cell).truth.stim_locations.(this_group).effect(:,this_loc_ID)*this_loc_power;
            end
        end
    end
    
    mu_bg = 1/experiment_setup.patched_neuron.background_rate;
    
    


    for i_trial = 1:number_of_trials
        
experiment_query_this_group.trials(i_trial).event_times=[];
experiment_query_this_group.trials(i_trial).assignments=[];
experiment_query_this_group.trials(i_trial).voltage_clamp=[];

        for i_cell = 1:number_of_cells
            k=stimuli_size(i_trial,i_cell);
            params_sim=experiment_setup.neurons(i_cell).truth;
            stim_threshold=experiment_setup.prior_info.induced_intensity.minimum_stim_threshold/experiment_setup.neurons(i_cell).truth.gain;
            if k > stim_threshold
                
                stim = experiment_setup.prior_info.current_template*k;
                [~, spikes]  =lif_glm_sim(stim,params_sim);
                %         plot(V_vect)
                if sum(spikes)>0
                    switch experiment_setup.prior_info.delay.type
                        case'gamma'
                        shape=(experiment_setup.prior_info.delay.mean^2)/...
                            (experiment_setup.prior_info.delay.std^2);
                        scale = experiment_setup.prior_info.delay.mean/shape;
                        delay_vec=round(gamrnd(shape,scale,[sum(spikes) 1]));
                    end
                    spikes_delay =find(spikes)+delay_vec';
                    for i = 1:sum(spikes)
                        if rand(1) < params_sim.PR
                            % censoring at time max:
                            if spikes_delay(i)<time_max
                                experiment_query_this_group.trials(i_trial).event_times=[experiment_query_this_group.trials(i_trial).event_times spikes_delay(i)];
                                experiment_query_this_group.trials(i_trial).assignments=[experiment_query_this_group.trials(i_trial).assignments i_cell];
                                
                            end
                        end
                        
                    end
                end
            end
        end
        % add background event:
        R = exprnd(mu_bg);
        while R < time_max
             experiment_query_this_group.trials(i_trial).event_times=[ experiment_query_this_group.trials(i_trial).event_times max(1,round(R))];
            experiment_query_this_group.trials(i_trial).assignments=[experiment_query_this_group.trials(i_trial).assignments 0];
            R = R+exprnd(mu_bg);
        end
        
    end
    % fprintf('%d trials', i_trial);
    % fprintf(' simulated;\n');

end

