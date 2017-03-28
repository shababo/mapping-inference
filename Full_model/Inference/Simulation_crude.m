%% Loading functions and Data generation

%%
     
    % Count the events in each amplitude bins:
    related_mpp_n=struct();
    for l = 1:n_trial
        if size(mpp_new(l).event_times,2) > 0
            indices = 1:size(mpp_new(l).event_times,2);
            related_mpp_n(l).amplitudes = mpp_new(l).amplitudes(indices);
            related_mpp_n(l).event_times = mpp_new(l).event_times(indices);
        else
            related_mpp_n(l).amplitudes = [];
            related_mpp_n(l).event_times = [];
        end
    end
    
	threshold=  quantile([related_mpp_n.amplitudes], (1/num_threshold)*[0:num_threshold]);
    
	for l = 1:n_trial
        for j = 1:num_threshold
            Y_g(l,j) = sum(related_mpp_n(l).amplitudes>threshold(j) & related_mpp_n(l).amplitudes<(threshold(j+1)+0.01));
        end
        
    end
    %%
    %---------------------------------------------------------%
    % Part III:
    % Fit the VB model to update the parameters:
    % Initialization:
    
    % Initialize starting values
    output_warm= struct([]);
    for j = 1:num_threshold
        output_warm(j).alpha = .1*ones(n_cell_local+1,1);
        output_warm(j).mu = zeros(n_cell_local+1,1);
        output_warm(j).s_sq = ones(n_cell_local+1,1);
    end
	
	
    X_temp =  stimuli_size_local;
    X_temp = [ones(size(X_temp,1),1) X_temp];
    for j = 1:size(Y_g,2)
        %fprintf('%d',j);
        Y_n = Y_g(:,j);
        if  sum(Y_n)==0 % prevent zero standard deviation
            hyperparam_sigma_n = 1;
        else
            hyperparam_sigma_n = std(Y_n);
        end
        
        hyperparam_p_connected = output_warm(j).alpha;
        hyperparam_eta =  output_warm(j).mu;
        hyperparam_sigma_s = sqrt(output_warm(j).s_sq);
        
        options = struct();
        options.alpha =  output_warm(j).alpha;
        options.mu= output_warm(j).mu;
        options.verbose= false;
        options.center= 0;
        
        [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(X_temp, Y_n, ...
            hyperparam_sigma_n, hyperparam_sigma_s, hyperparam_p_connected, ...
            hyperparam_eta, options);
        output(j).alpha = alpha_tmp;
        output(j).mu = mu_tmp;
        output(j).s_sq = s_sq_tmp;
    end
 
%% Evaluation
overall_connectivity = zeros(n_cell_local,1);
    overall_connectivity_w = zeros(n_cell_local,1);
    
    overall_mark = zeros(n_cell_local,1);
    normalized_constants = zeros(n_cell_local,1);
    for j = 1:num_threshold
        overall_connectivity = max(overall_connectivity, ...
     output(j).alpha(2:end));
        for i_cell = 1:n_cell_local
            if overall_connectivity(i_cell) == output(j).alpha(i_cell+1)
                overall_mark(i_cell)= (threshold(j)+threshold(j+1))/2;
            end
        end
     
    end
    
	
    %%
	%flnm=strcat('../../Data/Full_sim/Case', num2str(casenumber),'Crude.mat');
 
    %save(flnm, 'output','overall_mark','overall_connectivity');
    
