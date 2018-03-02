function [parameter_history] = fit_VI_dev3(...
    stim_size, mpp, background_rate, ...
    variational_params,prior_params,...
    inference_params,prior_info,spike_curves,varargin)
% Fit the delay for each cell:
% %%
% stim_size=designs_remained;
% background_rate=bg_rate;
% mpp=mpp_remained;
% spike_indicator=false;
%%
% et=[mpp.event_times];
% figure(1)
% histogram(et)
% xlim([0 300])
%%

if ~isempty(varargin) && ~isempty(varargin{1})
   spike_indicator= varargin{1};
else
    spike_indicator= false;
end
%%
intensity_grid=prior_info.induced_intensity.intensity_grid;
stim_scale=prior_info.induced_intensity.stim_scale;
minimum_stim_threshold=prior_info.induced_intensity.minimum_stim_threshold;

gamma_bound=struct;
gamma_bound.up=inference_params.bounds.PR(2);
gamma_bound.low=inference_params.bounds.PR(1);


gain_bound=struct;
gain_bound.up=inference_params.bounds.gain(2);
gain_bound.low=inference_params.bounds.gain(1);
epsilon=inference_params.convergence_threshold;
S=inference_params.MCsamples_for_gradient;
eta=inference_params.step_size;
eta_max=inference_params.step_size_max;
maxit=inference_params.maxit;
lklh_func=inference_params.likelihood;

delay_mu_bound=struct;
delay_mu_bound.up=60;
delay_mu_bound.low=0;

delay_sigma_bound=struct;
delay_sigma_bound.up=10;
delay_sigma_bound.low=0.1;


%stim_size=designs_remained; mpp= mpp_remained;

   
n_cell=size(stim_size,2);n_trial=size(stim_size,1);
n_grid=size(intensity_grid,2);

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


parameter_history=variational_params;
parameter_current=variational_params;

% find relevant trials for each cell 
relevant_trials = cell(n_cell,1);
for i_cell = 1:n_cell
    relevant_trials{i_cell}=find(stim_size(:,i_cell)>(minimum_stim_threshold/gain_bound.up));
end
changes=1;iteration = 1;
%%
while (change_history(iteration) > epsilon && iteration<maxit)
    %%
    parameter_history(iteration,:)=parameter_current;
    iteration=iteration+1;
    
    % precalculate some quantities 
    [gamma_constants] = get_logitnormal_constants([parameter_current(:).alpha],[parameter_current(:).beta]);
    [gain_constants] = get_logitnormal_constants([parameter_current(:).alpha_gain],[parameter_current(:).beta_gain]);
    [delay_mu_constants]=get_logitnormal_constants([parameter_current(:).alpha_m],[parameter_current(:).beta_m]);
    [delay_sigma_constants]=get_logitnormal_constants([parameter_current(:).alpha_s],[parameter_current(:).beta_s]);
    
    for s= 1:S
        % drawing samples 
        [logit_gamma, gamma_sample] = draw_spiked_logitnormal(...
            [parameter_current(:).alpha],[parameter_current(:).beta],gamma_bound,[parameter_current(:).p_logit],spike_indicator);
        [logit_gain, gain_sample] = draw_spiked_logitnormal(...
            [parameter_current(:).alpha_gain],[parameter_current(:).beta_gain],gain_bound);
        
        [logit_delay_mu, delay_mu_sample] = draw_spiked_logitnormal(...
            [parameter_current(:).alpha_m],[parameter_current(:).beta_m],delay_mu_bound);
        [logit_delay_sigma, delay_sigma_sample] = draw_spiked_logitnormal(...
            [parameter_current(:).alpha_s],[parameter_current(:).beta_s],delay_sigma_bound);
        
        
        delay_mu_sample_mat(:,s)=delay_mu_sample;
        delay_sigma_sample_mat(:,s)=delay_sigma_sample;
        gamma_sample_mat(:,s)=gamma_sample;
        gain_sample_mat(:,s)=gain_sample;
        
%         % Calculate log pdf of prior and variational dist.
        logprior(:,s)=calculate_logpdf_spiked_logitnormal_dev(prior_params,...
            logit_gamma,gamma_sample,logit_gain,gain_sample,...
            logit_delay_mu,delay_mu_sample,logit_delay_sigma,delay_sigma_sample,...
             delay_mu_bound,delay_sigma_bound,gamma_bound,gain_bound,spike_indicator);
       

        logvariational(:,s)=calculate_logpdf_spiked_logitnormal_dev(parameter_current,...
            logit_gamma,gamma_sample,logit_gain,gain_sample,...
            logit_delay_mu,delay_mu_sample,logit_delay_sigma,delay_sigma_sample,...
             delay_mu_bound,delay_sigma_bound,gamma_bound,gain_bound);
        
        
        [dqdp_logit(:,s),dqdalpha(:,s),dqdbeta(:,s)] = get_dvariational_dist(gamma_sample,logit_gamma,...
            [parameter_current(:).alpha],[parameter_current(:).beta],gamma_constants,[parameter_current(:).p_logit],spike_indicator);
        [~,dqdalpha_gain(:,s),dqdbeta_gain(:,s)] = get_dvariational_dist(gain_sample,logit_gain,...
            [parameter_current(:).alpha_gain],[parameter_current(:).beta_gain],gain_constants);
        [~,dqdalpha_s(:,s),dqdbeta_s(:,s)] = get_dvariational_dist(delay_sigma_sample,logit_delay_sigma,...
            [parameter_current(:).alpha_s],[parameter_current(:).beta_s],delay_sigma_constants);
        [~,dqdalpha_m(:,s),dqdbeta_m(:,s)] = get_dvariational_dist(delay_mu_sample,logit_delay_mu,...
            [parameter_current(:).alpha_m],[parameter_current(:).beta_m],delay_mu_constants);
        
    end
    
    %%
%     [mean(gamma_sample_mat)  [parameter_current(:).alpha]]
%     [var(gamma_sample_mat)  exp(2*[parameter_current(:).beta])]
%     [mean(gain_sample_mat)  [parameter_current(:).alpha_gain]]
%     [var(gain_sample_mat)  exp(2*[parameter_current(:).beta_gain])]
    

        % Only thing that needs to change if we use the new way to get
        % intensity: 
    [loglklh] = update_likelihood_dev3(gamma_sample_mat,gain_sample_mat,stim_size,mpp,...
            delay_mu_sample_mat,delay_sigma_sample_mat,...
            intensity_grid,minimum_stim_threshold,stim_scale,background_rate,relevant_trials,lklh_func,...
        spike_curves);
    
% scatter(gain_sample_mat,loglklh)
    [parameter_current, change_history(iteration)]= update_parameters_logitnormal_dev(...
        parameter_current,loglklh, logprior, logvariational,...
        dqdp_logit, dqdalpha, dqdbeta, dqdalpha_gain, dqdbeta_gain,...
         dqdalpha_m, dqdbeta_m, dqdalpha_s, dqdbeta_s,...
        iteration,eta,eta_max,spike_indicator);
% parameter_current

fprintf('Iteration %d; gain (alpha) %d; change %d; \n',iteration,parameter_current.alpha_gain, change_history(iteration))
end
fprintf('VI fitted after %d iteration;\n',iteration)
           
