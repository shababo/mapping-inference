%% ROI detection: Stage I
%  Directly Stimulating cell somas 

function output = ROI_VB_3D_grid(depth_data,stim_indices,depths,raw_traces)

%% build mpp

output= struct([]);

max_x = max(max(stim_indices(:,1,:)));
max_y = max(max(stim_indices(:,2,:)));
cell_coords = zeros(max_x*max_y,3);

    
% for map_i = 1:length(depth_data)
    
%     map_i
    
    mpp = depth_data;
    num_trials = length(mpp);

    %MAKE cell_coords
    unique_depths = unique(depths);
    for i = 1:max_x
        for j = 1:max_y
            for k = 1:length(unique_depths)
                cell_coords(sub2ind([max_x,max_y,length(unique_depths)],i,j,k),:)...
                    = [(i-1)*20-160 (j-1)*20-160 unique_depths(k)];
            end
        end
    end
    
    num_targets = size(stim_indices,1);
    
    trial_locations = zeros(num_trials,num_targets);
    for i = 1:num_trials
            
%         i
        trial_locations(i,:) = sub2ind([max_x,max_y,length(unique_depths)],stim_indices(:,1,i),stim_indices(:,2,i),ones(num_targets,1)*find(unique_depths == depths(i)));

    end
    %% params
    evoked_params.stim_start = 100;
    params.N = num_trials;



    %% Calculate summary statistics 

    % Use only the amp in the related regions 
    related_mpp = mpp;
    unrelated_mpp = mpp;

    % only consider a small time window for events of interest
    for i = 1:size(trial_locations,1)
        if size(mpp(i).times,2) > 0
            indices = mpp(i).times>evoked_params.stim_start & mpp(i).amp >= 0;
            related_mpp(i).amp = mpp(i).amp(indices);
            related_mpp(i).times = mpp(i).times(indices);
            unrelated_mpp(i).amp = mpp(i).amp(~indices);
            unrelated_mpp(i).times = mpp(i).times(~indices);
        end 
    end

    covariates = zeros(size(trial_locations,1), size(cell_coords,1));
    for i = 1:params.N
        covariates(i, trial_locations(i,:)) = 1;    
    end

    % With intercept, it is rank deficient. 
%     covariates_intercept = [ones(params.N,1) covariates];
    % rank(covariates)
    % size(covariates)


    %% Dividing the events by their amp
    % Now divide the events by quantiles of the amp 
    % We use overlapped regions to avoid separation due to 
    num_threshold=2;
    amplitude_threshold = quantile([related_mpp.amp], (1/num_threshold)*[0:num_threshold]);
    amp_related_count_trials = ones(size(trial_locations,1),num_threshold-1);
    for j = 1:(num_threshold-1)
        for i = 1:size(amp_related_count_trials,1)
            amp_related_count_trials(i,j) = sum(related_mpp(i).amp>amplitude_threshold(j) & related_mpp(i).amp<(amplitude_threshold(j+2)+0.01));
        end
    end

    %% The VB inference

    % params.A = A
    % params = struct;
    % params.A=A;
    params.coords=cell_coords;
    params.K = size(cell_coords,1); % num targets/cells
    % params.N=N; % num trials

    data=struct;
    data.stims = trial_locations;

    % Unknows: 
    params.eta = 1000 * ones(params.K,1);
    params.sigma_s = 1000*ones(params.K,1);
    if isempty(raw_traces)
        params.eta = ones(params.K,1);
        params.sigma_s = ones(params.K,1);
    end
%     if cell_type(map_i)
%         params.sigma_n = sqrt(2);
        params.A = [500 0 0 ; 0 500 0; 0 0 750];
        hyperparam_p_connected = .05*ones(params.K,1);
%     else
% %         params.sigma_n = sqrt(.75);
%         params.A = [100 0 0 ; 0 100 0; 0 0 100]*1;
%         hyperparam_p_connected = .01*ones(params.K,1);
%     end

%     params.t = 1:1:data_params.T;
    params.tau = 10;
    params.g = 1;
%     alpha_sum = sum(alpha_synapse(params.t,0,params.tau,-params.g));


    pi_kr = exp(-0.5*squareform(pdist(params.coords,'mahalanobis',params.A)).^2);

    pi_nk = zeros(params.N,params.K);
    for i = 1:params.N
        pi_nk(i,:) = min(1,sum(pi_kr(:,data.stims(i,:)),2)');
    end

    assignin('base','pi_nk',pi_nk)

    for j = 1:size(amp_related_count_trials,2)

        Y_n = amp_related_count_trials(:,j);
        if ~isempty(raw_traces)
            raw_traces_cent = bsxfun(@minus,raw_traces,median(raw_traces,2));
            Y_n = sum(raw_traces_cent(:,100:800),2);
            size(Y_n)
        end
        resp_sorted = sort(Y_n);%amp_related_count_trials(:,j);
        hyperparam_sigma_n = std(resp_sorted(1:ceil(length(resp_sorted*.80))));%sqrt(length(params.t))*params.sigma_n/abs(alpha_sum);
        assignin('base','Y_n',Y_n)


        alphas = zeros(params.K,1); %ones(params.K, 1) * alpha_0;
        mu = zeros(params.K, 1);
        s_sq = zeros(params.K,1);
        n_varbvs_samples = 10;
            options= struct();
        options.verbose= false;
        options.center= 1;

        for sample = 1:n_varbvs_samples
            [alpha_tmp, mu_tmp, s_sq_tmp] = run_varbvs_general(pi_nk>rand(params.N,params.K), Y_n, hyperparam_sigma_n, params.sigma_s(1), hyperparam_p_connected(1), params.eta,options);
            alphas = alphas+alpha_tmp/n_varbvs_samples;
            mu = mu+mu_tmp/n_varbvs_samples;
            s_sq = s_sq+s_sq_tmp/n_varbvs_samples;
        end
        output = struct();
        output.alpha = alphas;
        output.mu = mu;
        output.s_sq = s_sq;

        output.w_estimate = alphas.*mu;
    end
% end
%------------------------End of first stage-------------------------------------%

%% Visualize the Bayes estimates:

% figure(80)
% 
% for i = 1:num_layers
%     connected_neurons_ind = find(neuron_features(i).amplitude);
%     temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
%         -neuron_locations{i}(connected_neurons_ind,2),...
%         neuron_features(i).amplitude(connected_neurons_ind)*25);
%     set(temp,'MarkerFaceColor','k');
%    alpha(temp,0.8);
%     hold on
% end
% set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
% 
% %selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
% %cent = cell(size(amp_related_count_trials,2),1);
% 
%     xlim([20,460]);
%     ylim([-900,-400]);
% %     
% %      potential_neuron_grid = scatter(cell_coords(:,1),...
% %      -cell_coords(:,2),20,colormap(2,:),...
% %     'filled','d');
% %     set(potential_neuron_grid,'MarkerFaceColor','k');
% %     alpha(potential_neuron_grid,0.2);
% 
% for j = 1:size(amp_related_count_trials,2)
%     coef = output(map_i).ch1(j).alpha;
%     coef_thres = quantile(coef,0.98);
%     potential_neuron_grid = scatter(cell_coords(coef>coef_thres,1), -cell_coords(coef>coef_thres,2), amplitude_threshold(j+1)*25,'filled','o');
%     set(potential_neuron_grid,'MarkerFaceColor','r');
%     alpha(potential_neuron_grid,0.4);
%     hold on
%     
% end
% 
% hold off
% view(2)
% 
% %% 
% % Visualize the estiamted weights 
% 
% % merging the estimated weights
% w_estimate_merge = [];
% mu_merge = [];
% alpha_merge = [];
% for j = 1:size(amp_related_count_trials,2)
%     w_estimate_merge = [w_estimate_merge output(map_i).ch1(j).w_estimate];
%     mu_merge = [mu_merge output(map_i).ch1(j).mu];
%     alpha_merge = [alpha_merge output(map_i).ch1(j).alpha];
% end
% 
% %% Obtain true amp of the related neurons
% Z_amplitudes = all_amplitudes(neuron_in_region==1);
% 
% 
% 
% %% 
%    figure(11)
%    % find some cells that have non-zero amp, and some have zero
%    % amp
%    one_seq = 1:100;
%    nonzeros_seq = one_seq(Z_amplitudes(one_seq)>0);
%   chosen_ones = [1:5 nonzeros_seq(1:5)];
%    scale_y=0.5;
% for j = 1:length(chosen_ones)
%        
%     if Z_amplitudes(chosen_ones(j)) > 0
%         temp = line( amplitude_threshold(2:20),...
%             j/scale_y+2*w_estimate_merge(chosen_ones(j),:)/sum(w_estimate_merge(chosen_ones(j),:)),...
%             'LineWidth',3,...
%             'Color','k');
%         set(temp,'MarkerFaceColor','k');
%         %alpha(temp,0.8);
%         
%         %set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
%         
%         %selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
%         %cent = cell(size(amp_related_count_trials,2),1);
%         
%         x = [amplitude_threshold(1):.1:amplitude_threshold(20)];
%         norm = normpdf(x,Z_amplitudes(chosen_ones(j)),evoked_params.sigma_a);
%         
%         line(x,...
%             j/scale_y+norm,...
%             'LineWidth',4,...
%             'Color','r');
%         xlim([amplitude_threshold(1),amplitude_threshold(20)]);
%         ylim([-0.1,length(chosen_ones)/scale_y+1.1]);
%         %
%     else 
%         if sum(w_estimate_merge(chosen_ones(j),:))/2<0.1
%             scale = 1;
%         else 
%             scale = sum(w_estimate_merge(chosen_ones(j),:));
%         end
%         
%         temp = line( amplitude_threshold(2:20),...
%             j/scale_y+2*w_estimate_merge(chosen_ones(j),:)/scale,...
%             'LineWidth',3,...
%             'Color','b');
%         %set(temp,'MarkerFaceColor','b');
%         %alpha(temp, 0.3);
%     end
% end
% hold off

%% Waste code 


% %% Fit regression by amp 
% lmCount_related_amp=cell(size(amp_related_count_trials,2),1);
% %lmCount_related_amp_Robust=cell(size(amp_related_count_trials,2),1);
% for j = 1:size(amp_related_count_trials,2)
%     mdl_j=fitlm(covariates,amp_related_count_trials(:,j),'Intercept',false);
%     %EstCov = hac(mdl_j);
%     lmCount_related_amp{j}=mdl_j;
%     %lmCount_related_amp_Robust{j}=EstCov;
% end
% %% Visualize the estimates 
% % for j = 1:size(amp_related_count_trials,2)
% %     figure(j*17)
% %     histogram(lmCount_related_amp{j}.Coefficients.Estimate,30)
% %     
% %     % Check the quantiles and p-values 
% %     quantile(lmCount_related_amp{j}.Coefficients.Estimate,0.9)
% %     
% %     min(lmCount_related_amp{j}.Coefficients.Estimate(lmCount_related_amp{j}.Coefficients.pValue <0.05))
% %     
% % end
% figure(10)
% R=2500;
% 
% colormap = jet(2);
% colormap(1,:)= [1 1 1];
% colormap(2,:)= [1 0 0];
% for i = 1:num_layers
%     connected_neurons_ind = find(neuron_features(i).amplitude);
%     temp = scatter(neuron_locations{i}(connected_neurons_ind,1),...
%         -neuron_locations{i}(connected_neurons_ind,2),...
%         neuron_features(i).amplitude(connected_neurons_ind)*35);
%     set(temp,'MarkerFaceColor','k');
%     alpha(temp,0.8);
%     hold on
% end
% set(gca,'yticklabels',{'1200','1000','800','600','400','200','0'})
% 
% %selected_pixels = zeros(num_dense,num_dense,size(amp_related_count_trials,2));
% %cent = cell(size(amp_related_count_trials,2),1);
% 
%     xlim([20,460]);
%     ylim([-900,-400]);
% %     
% %      potential_neuron_grid = scatter(cell_coords(:,1),...
% %      -cell_coords(:,2),20,colormap(2,:),...
% %     'filled','d');
% %     set(potential_neuron_grid,'MarkerFaceColor','k');
% %     alpha(potential_neuron_grid,0.2);
% 
% for j = 1:size(amp_related_count_trials,2)
%     coef = lmCount_related_amp{j}.Coefficients.Estimate;
%     coef_thres = quantile(coef,0.98);
%     potential_neuron_grid = scatter(cell_coords(coef>coef_thres,1), -cell_coords(coef>coef_thres,2), amplitude_threshold(j+1)*35,'filled','o');
%     set(potential_neuron_grid,'MarkerFaceColor','r');
%     alpha(potential_neuron_grid,0.4);
% hold on
%    
% end
% 
% hold off
% %saveas(10,'../Data/Sites_to_stimulate.jpg')
% view(2)
