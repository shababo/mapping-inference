% Fit the gains based on the first spikes (soft assignments)
function [stats_conv, prob_first_spike_delayed] = fit_lifglm_v6(responses,stims,in_params,background_rate,v_reset_known,v_th_known,first_only,...
    loss_type,lower_bound,upper_bound, initial_values,delay_params,single_connection)
% responses is an N x T binary matrix of soft assignments where we have N trials
% stims is an N x T matrix of the stimulus timeseries which is scaled by power - if using a shape to
% fit the lif-glm then these should be scaled by the shape as well.
% in_params is struct with the field g which is 1/tau for this neuron

params.invlink = @invlink_sig;
params.dlink = @derlink_sig;
params.link = @link_sig;
params.dinvlink = @derinvlink_sig;
linkfunc = {params.link, params.dlink, params.invlink,params.dinvlink};

if first_only
    n_cell = size(stims,1);
    n_trial = size(stims,2);
    n_grid =size(stims,3);
    v_trace=zeros(n_cell, n_trial,n_grid);
    for i_cell = 1:n_cell
        for i_trial = 1:n_trial
            v_trace(i_cell,i_trial, 1) = 0;
            for i_t = 2:n_grid
                temp1=reshape(stims(i_cell, i_trial,1:(i_t-1)), [i_t-1,1]);
                temp2=reshape(exp( ((1:(i_t-1)) -i_t+1)*in_params.g(i_cell)), [1,i_t-1]);
                v_trace(i_cell, i_trial, i_t) = temp2*temp1;
            end
        end
    end
    
 if single_connection==true   
    betahat= zeros(size(initial_values,1)*(size(initial_values,2)/2),size(initial_values,2));
    loss_seq = zeros(size(betahat,1),1);
    ub =[upper_bound+1e-5 1+1e-5];
    lb = [lower_bound-1e-5 0];
    for i_cell = 1:(size(initial_values,2)/2)
        v_trace_single = v_trace(i_cell,:,:);
        
%         %% Debug:
%         delay_params.mean = 50;
%         xrange=gain_lower_bound + (1:20)/200;
% 
%         lklh=zeros(length(xrange),1);
%         
%         for i = 1:length(xrange)
%             lklh(i)=lif_glm_firstspike_loglikelihood_multi([xrange(i) 1], responses, ...
%                     v_trace_single,background_rate,v_th_known,linkfunc,loss_type,delay_params);
%             fprintf('%d',i);
%         end
%         
%         %% Draw spike time v.s. max stim
%         spike_time = [];
%         max_stim = [];
%         
%         for i_trial = 1:n_trial
%             spk_t=find(responses(i_trial,:));
%             if length(spk_t)>0
%         spike_time = [spk_t(1) spike_time];
%         max_stim = [max(stims(i_cell,i_trial,:)) max_stim];
%             end
%         end
%         plot(max_stim,spike_time,'.')
%         %%
%         figure(1)
% plot(xrange,lklh)
% figure(2)
% plot(reshape(v_trace_single(1,240,:), [n_grid 1]))
%         
%
        for i = 1:size(initial_values,1)
         x0=[initial_values(i,1) 1];
         
        obj_fun = @(x) lif_glm_firstspike_loglikelihood_multi(x, responses, v_trace_single,...
            background_rate,v_th_known,linkfunc,loss_type,delay_params);
        
        options = optimoptions('fmincon','Algorithm','interior-point',...
            'Display','iter',...
            'SpecifyObjectiveGradient',false,...
            'Diagnostics','off',...
            'UseParallel',true);
        %'GradObj','off',...
        [betahat(i+(i_cell-1)*size(initial_values,1),(1:2) +(i_cell-1)*2), loss_seq(i+(i_cell-1)*size(initial_values,1))] = ...
            fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);
        
   
        
        
        end
    end
    [~, min_idx] = min(loss_seq);
    stats_conv = betahat(min_idx,:);
else 
    betahat=initial_values;
    loss_seq = zeros(size(initial_values,1),1);
    ub = repmat([upper_bound+1e-5 1], [1 n_cell]);
    lb = repmat([lower_bound-1e-5 0], [1 n_cell]);
    
    for i=1:size(initial_values,1)
        x0=initial_values(i,:);
        obj_fun = @(x) lif_glm_firstspike_loglikelihood_multi(x, responses, v_trace,...
            background_rate,v_th_known,linkfunc,loss_type,delay_params);
        options = optimoptions('fmincon','Algorithm','interior-point',...
            'Display','iter',...
            'SpecifyObjectiveGradient',false,...
            'Diagnostics','off',...
            'UseParallel',true);
        %'GradObj','off',...
        
        
        [betahat(i,:), loss_seq(i)] = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],options);
    end
    [~, min_idx] = min(loss_seq);
    stats_conv = betahat(min_idx,:);
end
    
    
    % Prepare the probability:
    x= betahat(min_idx,:);
    prob_firing=zeros(n_cell,n_trial,n_grid);
    for i_cell = 1:n_cell
        for i_trial =1:n_trial
            for i_grid = 1:n_grid
                prob_firing(i_cell,i_trial,i_grid)=linkfunc{3}(x(2*(i_cell-1) +1)*v_trace(i_cell,i_trial,i_grid)-v_th_known(i_cell));
            end
        end
    end
    
    prob_first_spike=zeros(n_cell,n_trial,n_grid);
    for i_cell = 1:n_cell
        for i_trial =1:n_trial
            not_spike_prob=1;
            prob_first_spike(i_cell,i_trial,1)=prob_firing(i_cell,i_trial,1);
            for i_grid = 2:n_grid
                not_spike_prob = not_spike_prob*(1-prob_firing(i_cell,i_trial,i_grid -1));
                prob_first_spike(i_cell,i_trial, i_grid) =not_spike_prob*prob_firing(i_cell,i_trial,i_grid );
            end
        end
    end
    % Convoluted with delays:
    % Calculate the delay density on a grid
    
    if delay_params.delayed == 0
        prob_first_spike_delayed = prob_first_spike;% do nothing
    else
        delay_prob = zeros(delay_params.n_grid+1,1);
          if delay_params.type == 1 % normal
        delay_prob = normpdf((0:delay_params.n_grid),delay_params.mean,...
            delay_params.std);
    elseif delay_params.type == 2 %Gamma
        %shape 
        shape=(delay_params.mean^2)/(delay_params.std^2);
        %scale 
        scale = delay_params.mean/shape;
        delay_prob = gampdf((0:delay_params.n_grid),shape,scale);
    end
        % we approximate the probability with densities
        delay_prob = delay_prob/sum(delay_prob);
        min_delay = 0;
        max_delay = delay_params.n_grid;
        prob_first_spike_delayed = prob_first_spike;
        
        for i_cell = 1:n_cell
            for i_trial = 1:n_trial
                for i_grid = 1:n_grid
                    
                    idx_time = max(i_grid-max_delay,1): min(i_grid-min_delay,delay_params.n_grid);
                    idx_delay = -( (min(idx_time)-i_grid) : (max(idx_time)-i_grid))+1;
                    %M_grid_intensity{i_stimuli}(i_t)=delay_prob(idx_delay)*Intensity_grid{i_stimuli}(idx_time);
                    temp=0;
                    for i = 1:length(idx_time)
                        temp=temp+...
                            prob_first_spike(i_cell,i_trial,idx_time(i))*delay_prob(idx_delay(i));
                    end
                    prob_first_spike_delayed(i_cell,i_trial, i_grid) =temp;
                    
                end
            end
        end
    end


end