%% Simulate a two-neuron system with delays:
addpath(genpath('../../mapping-inference'));
bg_rate=1e-4;
%% Parameter range
location_range = 5:5:40;
gain_rg=[0.03 0.03];
PR_rg=[0.8 0.8];
dm_rg=[40  40];
dv_rg=[15 15];
rng(arg1)
%% Use real data (z-axis and xy )
path='../Environments/z_data.mat';
load(path);
axis_list={'x' 'y' 'z'};
result_z(1)=result_z(2);
path='../Environments/new_xy_data.mat';
load(path);
good_cell_list= [1:3 5:length(result_xy)];
for i_cell = good_cell_list
    
    for i_axis = 1:2
        ax=axis_list{i_axis};
        other_idx = 3-i_axis; % 2 or 1
        [unique_grid,ia,ic]=unique(abs(result_xy(i_cell).current_targ_pos(:,other_idx)));
        a_counts = accumarray(ic,1);
        [~,most_freq_loc]=max(a_counts);
        loc=unique_grid(most_freq_loc);
        this_axis_trials = find(ic==most_freq_loc);
        
        
        result_xy(i_cell).(ax)=struct;
        result_xy(i_cell).(ax).trials = this_axis_trials;
        result_xy(i_cell).(ax).stim_grid = result_xy(i_cell).current_targ_pos(this_axis_trials,i_axis);
        result_xy(i_cell).(ax).current = result_xy(i_cell).max_curr(this_axis_trials);
        result_xy(i_cell).(ax).power = result_xy(i_cell).curr_targ_power(this_axis_trials);
        
    end
    
end
result_xy(4)=result_xy(1); % replace the bad cell with some cell 1;

%% Start of the loops through 3 axes
%% Scale the current (devide by power and then by the max curr)
for i_ax = 3
    %%
    ax=axis_list{i_ax};
    if i_ax < 3
        
        n_cell = length(result_xy);
        clear('neurons')
        neurons(n_cell)=struct;
        for i_cell = 1:n_cell
            neurons(i_cell).stim_grid= result_xy(i_cell).(ax).stim_grid';
            %       neurons(i_cell).max_current=result_xy(i_cell).(ax).current';
            neurons(i_cell).power=result_xy(i_cell).(ax).power';
            
            neurons(i_cell).current=(result_xy(i_cell).(ax).current./result_xy(i_cell).(ax).power)';
            %         neurons(i_cell).scaled_current= neurons(i_cell).scaled_current/max( neurons(i_cell).scaled_current);
        end
        
    else
        n_cell = length(result_z);
        clear('neurons')
        neurons(n_cell)=struct;
        
        for i_cell = 1:n_cell
            neurons(i_cell).stim_grid=result_z(i_cell).current_z_pos;
            neurons(i_cell).current=result_z(i_cell).max_curr';
        end
        
    end
    %% Scale the current so the maximum is 1:
    % Calculate the measurement error (white noise)
    % Average the the recorded current on the same grid
    for i_cell = 1:n_cell
        [ugrid,~,ua]=unique(neurons(i_cell).stim_grid);
        neurons(i_cell).unique_grid = ugrid;
        neurons(i_cell).avg_current=zeros(length(ugrid),1);
        neurons(i_cell).sq_deviation=zeros(length(ugrid),1);
        
        for i_grid = 1:length(ugrid)
            indices = find(ua==i_grid);
            neurons(i_cell).avg_current(i_grid) =mean(neurons(i_cell).current(indices));
            neurons(i_cell).sq_deviation(i_grid) =var(neurons(i_cell).current(indices));
        end
        max_current = max(neurons(i_cell).avg_current);
        neurons(i_cell).scaled_current = neurons(i_cell).current/max_current;
        neurons(i_cell).scale=max_current;
    end
    
    %% Estimate the current noise:
    emp_current_sigma=zeros(n_cell,1);
    for i_cell = 1:n_cell
        emp_current_sigma(i_cell)=sqrt(mean([neurons(i_cell).sq_deviation])./[neurons(i_cell).scale]);
    end
    % plot(neurons(1).stim_grid,neurons(1).scaled_current)
    % 1 to ngrid +1
    found_shifts=zeros(n_cell,1);
    for i_cell = 1:n_cell
        [unique_stim,unique_ia,unique_ic]=unique(neurons(i_cell).stim_grid);
        mean_current =zeros(1,max(unique_ic));
        for i = 1:max(unique_ic)
            mean_current(i)= mean(neurons(i_cell).scaled_current(unique_ic==i));
        end
        
        n_grid = length( mean_current);
        [left]=prefixiso(n_grid,mean_current);
        [right]=prefixiso(n_grid,flip(mean_current));
        errorl=left.error(2:end);
        errorr=flip(right.error(2:end));
        totalerror=zeros(n_grid,1);
        for i=2:n_grid
            totalerror(i) = errorl(i-1)+errorr(i);
        end
        
        totalerror(1)=max(totalerror);
        [~,I_center]=min(totalerror);
        found_shifts(i_cell)=unique_stim(I_center);
    end
    %% Two-stage procedure
    neurons_original=neurons;
    for i_cell = 1:n_cell
        neurons(i_cell).stim_grid = neurons(i_cell).stim_grid-found_shifts(i_cell);
    end
    
    est_shape_func = fit([neurons(:).stim_grid]',[neurons(:).scaled_current]','smoothingspline','SmoothingParam',0.07);
    
    %% Generate a two-cell system:
    for i_loc = 1:length(location_range)
        loc=location_range(i_loc);
        neurons_2=struct;
        neurons_2(1)=struct;
        neurons_2(2)=struct;
        neurons_2(1).gain=0.06;
        % neurons_2(1).gain=0.05;
        
        neurons_2(1).PR=0;
        neurons_2(1).location=0;
        neurons_2(1).delay_mean=30;
        neurons_2(1).delay_var=60;
        
        neurons_2(2).gain=unifrnd(gain_rg(1),gain_rg(2));
        % neurons_2(2).gain=0.03;
        neurons_2(2).PR=unifrnd(PR_rg(1),PR_rg(2));
        neurons_2(2).location=loc;
        neurons_2(2).delay_mean=unifrnd(dm_rg(1),dm_rg(2));
        neurons_2(2).delay_var=unifrnd(dv_rg(1),dv_rg(2));
        
        loc_target= [neurons_2(:).location];
        
        %% Calculate the shape for each cell per stim locaiton
        % quantile_probs=[0.25 0.75];
        for i=1:2
            rel_loc = loc_target-neurons_2(i).location;
            neurons_2(i).z=struct;
            neurons_2(i).z.received=max(0,est_shape_func(rel_loc));
            %    neurons_2(i).z.received_low=max(0,z_quantile_function{1}(rel_loc));
            %    neurons_2(i).z.received_up=max(0,z_quantile_function{2}(rel_loc));
        end
        %% Design trials:
        Nrep=10;
        
        powers=[30 60 90];
        Ntrials = Nrep*length(loc_target)*length(powers);
        clear('mpp');
        mpp(Ntrials)=struct;
        for i_loc = 1:length(loc_target)
            for i_pow = 1:length(powers)
                for i_rep = 1:Nrep
                    n= (i_loc-1)*length(powers)*Nrep + (i_pow-1)*Nrep+i_rep;
                    mpp(n).power=powers(i_pow);
                    mpp(n).location=i_loc;
                    mpp(n).location_coord=loc_target(i_loc);
                    mpp(n).stimulation = mpp(n).power*[neurons_2(1).z.received(i_loc) neurons_2(2).z.received(i_loc)];
                end
            end
            
        end
        
        %% Simulate the point processes
        single_patch_path = '../Data/more_cells.mat';
        
        % calculate the respon se curve for simulation
        spike_curves=get_spike_curves(single_patch_path);
        
        bg_rate=1e-5;
        Tmax=300;
        for i= 1:Ntrials
            spike_times = [];
            event_times = [];
            assignments=[];
            eff_size=[];
            for i_neuron = 1:length(neurons_2)
                actual_stim= neurons_2(i_neuron).gain*mpp(i).stimulation(i_neuron);
                [~, Ia]=min(abs(actual_stim - spike_curves.current));
                spike_param = struct;
                spike_param.mean=spike_curves.mean(Ia);
                spike_param.sd=spike_curves.sd(Ia);
                spike_one = normrnd(spike_param.mean,spike_param.sd);
                delay_one = normrnd(neurons_2(i_neuron).delay_mean,sqrt(neurons_2(i_neuron).delay_var));
                % truncate the event if it is larger than Tmax
                if ((spike_one+delay_one)<Tmax) & (unifrnd(0,1)<neurons_2(i_neuron).PR)
                    spike_times = [spike_times spike_one];
                    assignments=[assignments i_neuron];
                    event_times = [event_times spike_one+delay_one];
                    eff_size=[eff_size mpp(i).stimulation(i_neuron)];
                end
            end
            % injecting background events:
            bg_prob=bg_rate*Tmax;
            bg_yes = binornd(1,bg_prob);
            if bg_yes
                spike_one=unifrnd(0,Tmax);
                spike_times = [spike_times spike_one];
                assignments=[assignments 0];
                event_times = [event_times spike_one];
                eff_size=[eff_size 0];
            end
            
            mpp(i).event_times = event_times;
            mpp(i).spike_times = spike_times;
            mpp(i).assignments = assignments;
            mpp(i).eff_size = eff_size;
        end
        % clear('spike_curves');
     
        %% Now apply the SVI
        
        
        disp('Test get_experiment_setup')
        load('./power-calibration.mat'); % somehow needs to load this manually
        experiment_setup=get_experiment_setup('szchen-sim-hab');
        %%
        group_profile=experiment_setup.groups.undefined;
        inference_params=group_profile.inference_params;
        % inference_params.MCsamples_for_gradient=100;
        % inference_params.convergence_threshold=1e-3;
        % inference_params.step_size_max=1;
        fitted_neurons=struct;
        
        neurons=struct;
        
        
        prior_info=experiment_setup.prior_info;
        temp_output=struct;
        % Non-informative initialization:
        temp_output.PR_params=struct;
        temp_output.PR_params.pi_logit=-Inf*ones(2,1);
        temp_output.PR_params.alpha=0*ones(2,1);
        temp_output.PR_params.beta=1*ones(2,1);
        
        temp_output.gain_params=temp_output.PR_params;
        temp_output.delay_mu_params=temp_output.PR_params;
        temp_output.delay_sigma_params=temp_output.PR_params;
        
        number_of_stim_cells=2;
        variational_params=struct;
        temp=num2cell(temp_output.PR_params.pi_logit); [variational_params(1:number_of_stim_cells).p_logit]=temp{:};
        temp=num2cell(temp_output.PR_params.alpha); [variational_params(1:number_of_stim_cells).alpha]=temp{:};
        temp=num2cell(temp_output.PR_params.beta); [variational_params(1:number_of_stim_cells).beta]=temp{:};
        temp=num2cell(temp_output.gain_params.alpha); [variational_params(1:number_of_stim_cells).alpha_gain]=temp{:};
        temp=num2cell(temp_output.gain_params.beta); [variational_params(1:number_of_stim_cells).beta_gain]=temp{:};
        temp=num2cell(temp_output.delay_mu_params.alpha); [variational_params(1:number_of_stim_cells).alpha_m]=temp{:};
        temp=num2cell(temp_output.delay_mu_params.beta); [variational_params(1:number_of_stim_cells).beta_m]=temp{:};
        temp=num2cell(temp_output.delay_sigma_params.alpha); [variational_params(1:number_of_stim_cells).alpha_s]=temp{:};
        temp=num2cell(temp_output.delay_sigma_params.beta); [variational_params(1:number_of_stim_cells).beta_s]=temp{:};
        
        % save('tmp.mat')
        prior_params=variational_params;
        
        %prior_params=variational_params_path(max(iter-num_trace_back,1),cell_list);
        
        bg_rate=1e-5; %experiment_setup.patched_neuron.background_rate
        
        %% Fit the model with the fit_VI
        stim_size=reshape([mpp(:).stimulation], [2 Ntrials])';
        
        [parameter_history,loglklh_rec] = fit_VI(...
            stim_size, mpp, bg_rate,...
            variational_params,prior_params,...
            inference_params,prior_info);
        
        %%
        for this_cell = 1:length(neurons_2)
            
            quantile_prob=group_profile.regroup_func_params.quantile_prob;
            current_params=reformat_to_neurons(parameter_history(end, this_cell),'PR','spiked_logit_normal');
            %     group_profile=experiment_setup.groups.(this_neighbourhood.neurons(i_cell).group_ID);
            bounds= group_profile.inference_params.bounds.PR;
            fitted_neurons(this_cell).PR_params=...
                calculate_posterior(current_params,bounds,quantile_prob);
            
            current_params=reformat_to_neurons(parameter_history(end, this_cell),'gain','spiked_logit_normal');
            bounds= group_profile.inference_params.bounds.gain;
            fitted_neurons(this_cell).gain_params=calculate_posterior(...
                current_params,bounds,quantile_prob);
            
            current_params=reformat_to_neurons(parameter_history(end, this_cell),'delay_mu','spiked_logit_normal');
            bounds= group_profile.inference_params.bounds.delay_mu;
            fitted_neurons(this_cell).delay_mu=calculate_posterior(...
                current_params,bounds,quantile_prob);
            
            current_params=reformat_to_neurons(parameter_history(end, this_cell),'delay_sigma','spiked_logit_normal');
            bounds=group_profile.inference_params.bounds.delay_sigma;
            fitted_neurons(this_cell).delay_sigma=calculate_posterior(...
                current_params,bounds,quantile_prob);
        end
        %% Reformat:
        clear('delay_params_svi')
        delay_params_svi(2)=struct;
        for i_neuron = 1:length(neurons_2)
            delay_params_svi(i_neuron).var=(fitted_neurons(i_neuron).delay_sigma.mean)^2;
            delay_params_svi(i_neuron).mean=(fitted_neurons(i_neuron).delay_mu.mean);
            delay_params_svi(i_neuron).gain=fitted_neurons(i_neuron).gain_params.mean;
            delay_params_svi(i_neuron).mpp=mpp;
            delay_params_svi(i_neuron).stim=stim_size(:,i_neuron);
            delay_params_svi(i_neuron).PR=fitted_neurons(i_neuron).PR_params.mean;
            
        end
        
        file_path = ['./Data/Sim/Axis' ax 'Loc' num2str(loc) 'Sim' num2str(arg1) '.mat'];
        save(file_path,'neurons_2', 'delay_params_svi','mpp')
    end
end