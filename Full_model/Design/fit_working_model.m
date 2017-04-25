function [output, Y_g,X_g] = fit_working_model(...
    i_batch,locations_this_batch,mpp_new,Y_g,X_g,output_old,num_threshold,evoked_params,length_memory,k_offset,...
    mark,pi_dense_local)

%total_events = [mpp_new.assignments];
%precells = local_index(local_connected);

% Count the events in each amplitude bins:
num_this_batch = size(mpp_new,2);
related_mpp_n=struct();
for l = 1:num_this_batch
    if size(mpp_new(l).event_times,2) > 0
        indices = mpp_new(l).event_times>evoked_params.stim_start & ...
            mpp_new(l).event_times< evoked_params.stim_end;
            %(20+evoked_params.stim_start);
        related_mpp_n(l).amplitudes = mpp_new(l).amplitudes(indices);
        related_mpp_n(l).event_times = mpp_new(l).event_times(indices)-evoked_params.stim_start;
    else
        related_mpp_n(l).amplitudes = [];
        related_mpp_n(l).event_times = [];
    end
end

if i_batch == 1
    quans =  (1/num_threshold)*[0:num_threshold];
    quans(end) = quans(end)-0.01; % prevent heavy-tail
    if mark == 0
        threshold=  quantile([related_mpp_n.amplitudes],quans);
        % Alternatively, use an evenly-spaced grid 
        %threshold= linspace(min([related_mpp_n.amplitudes]) ,max([related_mpp_n.amplitudes]),num_threshold+1);
    elseif mark == 1
        threshold=  quantile([related_mpp_n.event_times], quans);
    end
else
    threshold = output_old(1).threshold;
    num_threshold = size(output_old,2);
end
Y_new = zeros(num_this_batch,num_threshold);
for l = 1:num_this_batch
    for j = 1:num_threshold
        if mark == 0
            Y_new(l,j) = sum(related_mpp_n(l).amplitudes>threshold(j) & related_mpp_n(l).amplitudes<(threshold(j+1)+0.01));
        elseif mark==1
            Y_new(l,j) = sum(related_mpp_n(l).event_times >threshold(j) & related_mpp_n(l).event_times <(threshold(j+1)+0.01));
        end
    end
end

Y_g = [Y_g; Y_new];

X_next = zeros(size(pi_dense_local,1), num_this_batch);
 for l = 1:num_this_batch
        X_next(:,l) =  sum(pi_dense_local(:,locations_this_batch(l,:)),2);
 end
    
X_g =[X_g; X_next'];

%---------------------------------------------------------%
% Fit the VB model to update the parameters:
% Initialization:
output_warm = output_old;
if i_batch == 1
    ind_samples = 1:size(X_g,1);
else
    ind_samples = max(1,size(X_g,1)-length_memory):size(X_g,1);
end
 X_temp = k_offset.*X_g(ind_samples,:);
X_temp = [ones(size(X_temp,1),1) X_temp];

for j = 1:size(Y_g,2)
    %fprintf('%d',j);
   Y_n = Y_g(ind_samples,j);
    
   
    if  sum(Y_n)==0 % prevent zero standard deviation
        hyperparam_sigma_n = 1;
    else
        hyperparam_sigma_n = std(Y_n);
    end
    
    
    hyperparam_p_connected = output_warm(j).alpha; % No initial values for alpha
    hyperparam_eta =  output_warm(j).mu;
    hyperparam_sigma_s = sqrt(output_warm(j).s_sq);
    %ones(n_cell_local+1,1);
    %
    %
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
    output(j).threshold= threshold;
    
end
