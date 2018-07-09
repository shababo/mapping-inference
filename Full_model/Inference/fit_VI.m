function [parameter_history, elbo_rec,loglklh_rec] = fit_VI(...
    trials,neurons, background_rate, ...
    variational_params,prior_params,...
    inference_params,prior_info)
%%
n_cell=length(neurons);
n_trial = length(trials);

epsilon=inference_params.convergence_threshold;
S=inference_params.MCsamples_for_gradient;
eta=inference_params.step_size;
eta_max=inference_params.step_size_max;
maxit=inference_params.maxit;
lklh_func=inference_params.likelihood;

change_history=epsilon+1;

%%
boundary_params = prior_info.prior_parameters.boundary_params;
GP_params=prior_info.prior_parameters.GP_params;
spike_curves=prior_info.induced_intensity;
parameter_current=variational_params;
  

iteration = 1;loglklh_rec=[];
%%
% tic;
% time_rec=zeros(10,1);
tic;
while (change_history(iteration) > epsilon && iteration<maxit)
    %%
    tstart=toc;
    parameter_history(iteration,:)=parameter_current;
    iteration=iteration+1;
    adjusted_location=zeros(n_cell,3);
    loglklh=zeros(S,1);logprior=zeros(S,1);logvariational=zeros(S,1);
    log_prior_GP=zeros(n_cell,S);cell_shapes = cell([n_cell S]);shapes_indices= cell([n_cell S]);
    %% 
   
    
    %% Calculate the stim size using samples from GP and the shifts
    vsam=cell([S 1]);rsam=cell([S 1]);
    for s=1:S
        [variational_samples,raw_samples] = draw_samples_from_var_dist(parameter_current);
        vsam{s}=variational_samples;rsam{s}=raw_samples;
    end
    
%     tic;
%     ts=zeros(10,1);
    for s = 1:S
        
        variational_samples=vsam{s};
        raw_samples=rsam{s};
        
        for i_cell = 1:n_cell
%             tb=toc;
   
            shifts=[variational_samples(i_cell).shift_x variational_samples(i_cell).shift_y variational_samples(i_cell).shift_z];
            adjusted_location(i_cell,:)=neurons(i_cell).location-shifts;
%                       te=toc;
%     ts(1)=ts(1)+te-tb;
   
            relevant_trials=[];
            relevant_stim_locations=zeros(0,3);
            trial_and_stim=zeros(0,2);
             
            
            this_adjusted_loc=adjusted_location(i_cell,:);
            for i_trial = 1:length(trials)
                relevant_flag = false;
                for i_loc = 1:size(trials(i_trial).locations,1)
                    rel_position=trials(i_trial).locations(i_loc,:)-this_adjusted_loc;
                    this_relevant_flag =check_in_boundary(rel_position,boundary_params);
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
    
%     tb=toc;
            if ~isempty(relevant_trials)
                % Draw a GP on these positions:
                
                [GP_tmp] = draw_3D_GP(relevant_stim_locations,1, GP_params);
                cell_shapes{i_cell,s} = GP_tmp.full.samples; % n_loc by S vector
                shapes_indices{i_cell,s}=trial_and_stim;
                log_prior_GP(i_cell,s)=GP_tmp.full.loglklh;
                
            end
%              te=toc;
%     ts(2)=ts(2)+te-tb;
%     
        end
    end
      %
    for s= 1:S
        stim_size=zeros(n_trial,n_cell);
%         [variational_samples,raw_samples] = draw_samples_from_var_dist(parameter_current);
        % Store the stim size in stim all
           variational_samples=vsam{s}; raw_samples=rsam{s};
        for i_cell = 1:n_cell
            this_cell_shape=cell_shapes{i_cell,s};
            trial_and_stim=shapes_indices{i_cell,s};
            for i_rel = 1:length(this_cell_shape)
                i_trial = trial_and_stim(i_rel,1);
                i_loc = trial_and_stim(i_rel,2);
                if length(trials(i_trial).power_levels)==1
                    power_tmp =  trials(i_trial).power_levels;
                else
                    power_tmp = trials(i_trial).power_levels(i_loc);
                end
                stim_size(i_trial,i_cell)=stim_size(i_trial,i_cell)+ power_tmp*this_cell_shape(i_rel);
            end
        end
      
        
        %
        
        % Draw samples from the GP on this corrected grid
        logprior(s)=get_logdistribution(variational_samples,raw_samples,prior_params)+...
            sum(log_prior_GP(:,s));
        % NOTE: need to add the prior distribution of GP
        
        %         ts=toc;
        %         % Calculate log pdf of prior and variational dist.
        %         logprior(:,s)=calculate_logpdf_spiked_logitnormal(prior_params,...
        %             logit_gamma,gamma_sample,logit_gain,gain_sample,...
        %             logit_delay_mu,delay_mu_sample,logit_delay_sigma,delay_sigma_sample,...
        %              delay_mu_bound,delay_sigma_bound,gamma_bound,gain_bound,spike_indicator);
        
        
        logvariational(s)=get_logdistribution(variational_samples,raw_samples,parameter_current);
        
        %         logvariational(:,s)=calculate_logpdf_spiked_logitnormal(parameter_current,...
        %             logit_gamma,gamma_sample,logit_gain,gain_sample,...
        %             logit_delay_mu,delay_mu_sample,logit_delay_sigma,delay_sigma_sample,...
        %              delay_mu_bound,delay_sigma_bound,gamma_bound,gain_bound);
        %         te=toc;        time_rec(2)=time_rec(2)+te-ts;
        
        %         ts=toc;
        
        % Calcualte likelihood:
%           tb=toc;
      
        [loglklh(s)] = update_likelihood(stim_size,trials,...
            variational_samples, background_rate,lklh_func,spike_curves);
%          te=toc;
%     ts(3)=ts(3)+te-tb;
    
        lklhweight=logprior(s)+loglklh(s)-logvariational(s);
        
%         tb=toc;
        this_gradient=get_variational_gradient(variational_samples,raw_samples, parameter_current);
        this_gradient=get_variate_control(lklhweight,this_gradient);
%          te=toc;
%     ts(4)=ts(4)+te-tb;
        if exist('gradients')
            gradients(s,:) = this_gradient;
        else
            gradients(s,:)=this_gradient;
        end
        
    end
    %%
    new_gradient=sum_gradient(gradients,eta,eta_max,iteration);
    
    [parameter_current, change_history(iteration)]=incorporate_gradient(parameter_current, new_gradient);
    
    loglklh_rec(iteration)=mean(mean(loglklh));
    elbo_rec(iteration)=mean(logprior+loglklh-logvariational);
    tdiff=toc-tstart;
      fprintf('Iteration %d; change %d; ELBO %d; time %d; \n',iteration,change_history(iteration),elbo_rec(iteration),tdiff)
end
fprintf('VI fitted after %d iteration;\n',iteration)

