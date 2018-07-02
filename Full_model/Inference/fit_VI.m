function [parameter_history, loglklh_rec] = fit_VI(...
    trials,neurons, background_rate, ...
    variational_params,prior_params,...
    inference_params,prior_info,varargin)
%%

if ~isempty(varargin) && ~isempty(varargin{1})
    zero_delay_indicator= varargin{1};
else
    zero_delay_indicator= false;
end
% if ~isempty(varargin) && ~isempty(varargin{2})
%    spike_indicator= varargin{2};
% else
spike_indicator= false;
% end

if zero_delay_indicator
    variational_params.alpha_s=-200;
    variational_params.alpha_m=-200;
    
    prior_params.alpha_s=-200;
    prior_params.alpha_m=-200;
end

%%
% intensity_grid=prior_info.induced_intensity.intensity_grid;
% stim_scale=prior_info.induced_intensity.stim_scale;
spike_curves=prior_info.induced_intensity;
minimum_stim_threshold=spike_curves.minimum_stim_threshold;
gamma_bound=struct;
gamma_bound.up=inference_params.bounds.PR(2);
gamma_bound.low=inference_params.bounds.PR(1);

gain_bound=struct;
gain_bound.up=inference_params.bounds.gain(2);
gain_bound.low=inference_params.bounds.gain(1);

delay_mu_bound=struct;
delay_mu_bound.up=inference_params.bounds.delay_mu(2);
delay_mu_bound.low=inference_params.bounds.delay_mu(1);

delay_sigma_bound=struct;
delay_sigma_bound.up=inference_params.bounds.delay_sigma(2);
delay_sigma_bound.low=inference_params.bounds.delay_sigma(1);

epsilon=inference_params.convergence_threshold;
S=inference_params.MCsamples_for_gradient;
eta=inference_params.step_size;
eta_max=inference_params.step_size_max;
maxit=inference_params.maxit;
lklh_func=inference_params.likelihood;

n_cell=length(neurons);

Tmax=spike_curves.time_max;

% initialize storages
sum_of_logs=zeros(S,1);logvariational=zeros(n_cell,S);
logprior=zeros(n_cell,S);loglklh=zeros(n_cell,S);

dqdp_logit=zeros(n_cell,S);
dqdalpha=zeros(n_cell,S);dqdbeta=zeros(n_cell,S);
dqdalpha_gain=zeros(n_cell,S);dqdbeta_gain=zeros(n_cell,S);
dqdalpha_m=zeros(n_cell,S);dqdbeta_m=zeros(n_cell,S);
dqdalpha_s=zeros(n_cell,S);dqdbeta_s=zeros(n_cell,S);

change_history=epsilon+1;
gamma_sample_mat=zeros(n_cell,S);
gain_sample_mat=zeros(n_cell,S);
delay_mu_sample_mat=zeros(n_cell,S);
delay_sigma_sample_mat=zeros(n_cell,S);

%%
boundary_params = prior_info.boundary_params;
parameter_current=variational_params;

% % find relevant trials for each cell
% relevant_trials = cell(n_cell,1);
% for i_cell = 1:n_cell
%     relevant_trials{i_cell}=find(stim_size(:,i_cell)>(minimum_stim_threshold/gain_bound.up));
% end
iteration = 1;loglklh_rec=[];
%%
% tic;
% time_rec=zeros(10,1);
while (change_history(iteration) > epsilon && iteration<maxit)
    %%
    parameter_history(iteration,:)=parameter_current;
    iteration=iteration+1;
    log_prior_GP=zeros(n_cell,1);
    adjusted_location=zeros(n_cell,3);
    for s= 1:S
        % drawing samples
        [variational_samples,raw_samples] = draw_samples_from_var_dist(parameter_current);
        %% Calculate the stim size using samples from GP and the shifts
        for i_cell = 1:length(variational_samples) 
        shifts=[variational_samples(i_cell).shift_x variational_samples(i_cell).shift_y variational_samples(i_cell).shift_z];
        adjusted_location(i_cell,:)=neurons(i_cell).location-shifts;
        
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
                if relevant_flag
                relevant_trials = [relevant_trials i_trial];
                end
            end
        end
        
        % Draw a GP on these positions:
        [this_cell_shape] = draw_3D_GP(relevant_stim_locations,n_shapes, GP_params);
        
        % Store the stim size in stim all
        for i_rel = 1:size(relavant_stim_locations,1)
            i_trial = trial_and_stim(i_rel,1);
            i_loc = trial_and_stim(i_rel,2);
            stim_all(i_trial,i_cell)= trials(i_trial).power(i_loc)*this_cell_shape(i_rel);
        end
        log_prior_GP(i_cell)=this_cell_shape.loglklh;
        
        
        end
        
        %%
        
        % Draw samples from the GP on this corrected grid
        
        logprior(s)=get_logdistribution(variational_samples,prior_params)+...
            sum(log_prior_GP);
        % NOTE: need to add the prior distribution of GP
        
        %         ts=toc;
        %         % Calculate log pdf of prior and variational dist.
        %         logprior(:,s)=calculate_logpdf_spiked_logitnormal(prior_params,...
        %             logit_gamma,gamma_sample,logit_gain,gain_sample,...
        %             logit_delay_mu,delay_mu_sample,logit_delay_sigma,delay_sigma_sample,...
        %              delay_mu_bound,delay_sigma_bound,gamma_bound,gain_bound,spike_indicator);
        
        
        logvariational(s)=get_logdistribution(variational_samples,parameter_current);
        
        %         logvariational(:,s)=calculate_logpdf_spiked_logitnormal(parameter_current,...
        %             logit_gamma,gamma_sample,logit_gain,gain_sample,...
        %             logit_delay_mu,delay_mu_sample,logit_delay_sigma,delay_sigma_sample,...
        %              delay_mu_bound,delay_sigma_bound,gamma_bound,gain_bound);
        %         te=toc;        time_rec(2)=time_rec(2)+te-ts;
        
        %         ts=toc;
        
        % Calcualte likelihood:
        
        [loglklh(s)] = update_likelihood(stim_size,mpp,...
            gamma_sample_mat,gain_sample_mat,delay_mu_sample_mat,delay_sigma_sample_mat,...
            minimum_stim_threshold,background_rate,relevant_trials,lklh_func,...
            spike_curves,Tmax);
        
        
        lklhweight=logprior(s)+loglklh(s)-logvariational(s);
        
        this_gradient=get_variational_gradient(variational_samples,parameter_current);
        this_gradient=get_variate_control(lklhweight,this_gradient);
        
        if exist('gradients')
            gradients(s,:) = this_gradient;
        else
            gradients(s,:)=this_gradient;
        end
        
    end
    
    new_gradient=sum_gradient(gradients,eta,eta_max,iteration);
    
    [parameter_current, change_history(iteration)]=incorporate_gradient(parameter_current, new_gradient);
    
    loglklh_rec(iteration)=mean(mean(loglklh));
    elbo_rec(iteration)=mean(logprior+loglklh-logvariational);
    
    
    % Only thing that needs to change if we use the new way to get
    % intensity:
    %         ts=toc;
    [loglklh] = update_likelihood(stim_size,mpp,...
        gamma_sample_mat,gain_sample_mat,delay_mu_sample_mat,delay_sigma_sample_mat,...
        minimum_stim_threshold,background_rate,relevant_trials,lklh_func,...
        spike_curves,Tmax);
    %         te=toc;        time_rec(4)=time_rec(4)+te-ts;
    % scatter(gain_sample_mat,loglklh)
    % ts=toc;
    [parameter_current, change_history(iteration)]= update_parameters_logitnormal(...
        parameter_current,loglklh, logprior, logvariational,...
        dqdp_logit, dqdalpha, dqdbeta, dqdalpha_gain, dqdbeta_gain,...
        dqdalpha_m, dqdbeta_m, dqdalpha_s, dqdbeta_s,...
        iteration,eta,eta_max,spike_indicator);
    %         te=toc;        time_rec(5)=time_rec(5)+te-ts;
    
    % parameter_current
    loglklh_rec(iteration)=mean(mean(loglklh));
    % change_history(iteration)=abs(loglklh_rec(iteration)-loglklh_rec(iteration-1))/abs(loglklh_rec(iteration));
    % parameter_current.alpha_m=-5;
    % parameter_current.beta_m=-1;
    % parameter_current.beta_s=-1;
    % parameter_current.alpha_s=-1;
    fprintf('Iteration %d; change %d; \n',iteration,change_history(iteration))
end
fprintf('VI fitted after %d iteration;\n',iteration)

