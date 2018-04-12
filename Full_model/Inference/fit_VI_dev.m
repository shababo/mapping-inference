function [parameter_history, loglklh_rec] = fit_VI_dev(...
 mpp, background_rate,z_dictionary,...
            variational_params,prior_params,...
            inference_params,prior_info,varargin)
        
% Fit the delay for each cell:
% %%
% power_size=designs_remained;
% background_rate=bg_rate;
% mpp=mpp_remained;
% spike_indicator=false;
%%
%   background_rate=experiment_setup.patched_neuron.background_rate;
% inference_params=  group_profile.inference_params;
%  spike_indicator= false;
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

%% Calculate power size:
if length(mpp(1).power)==length(prior_params)
    power_size=reshape([mpp(:).power], [2 length(mpp)])';
else
    power_size=([mpp(:).power]')*ones(1,length(prior_params));
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

n_cell=size(power_size,2);n_trial=size(power_size,1);

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
shape_sample_mat=zeros(n_cell,S);

parameter_history=variational_params;
parameter_current=variational_params;

% find relevant trials for each cell 
relevant_trials = cell(n_cell,1);
for i_cell = 1:n_cell
    relevant_trials{i_cell}=find(power_size(:,i_cell)>(minimum_stim_threshold/gain_bound.up));
end
iteration = 1;loglklh_rec=[];
%%
% tic;
% time_rec=zeros(10,1);
while (change_history(iteration) > epsilon && iteration<maxit)
    %%
    parameter_history(iteration,:)=parameter_current;
%     parameter_current=parameter_history(end,:);
  
    iteration=iteration+1;
    
    % precalculate some quantities 
    [gamma_constants] = get_logitnormal_constants([parameter_current(:).alpha],[parameter_current(:).beta]);
    [gain_constants] = get_logitnormal_constants([parameter_current(:).alpha_gain],[parameter_current(:).beta_gain]);
    [delay_mu_constants]=get_logitnormal_constants([parameter_current(:).alpha_m],[parameter_current(:).beta_m]);
    [delay_sigma_constants]=get_logitnormal_constants([parameter_current(:).alpha_s],[parameter_current(:).beta_s]);
    
    for s= 1:S
        % drawing samples 
%         t1s=toc;
        [logit_gamma, gamma_sample] = draw_spiked_logitnormal(...
            [parameter_current(:).alpha],[parameter_current(:).beta],gamma_bound,[parameter_current(:).p_logit],spike_indicator);
        [logit_gain, gain_sample] = draw_spiked_logitnormal(...
            [parameter_current(:).alpha_gain],[parameter_current(:).beta_gain],gain_bound);
        
        [logit_delay_mu, delay_mu_sample] = draw_spiked_logitnormal(...
            [parameter_current(:).alpha_m],[parameter_current(:).beta_m],delay_mu_bound);
        [logit_delay_sigma, delay_sigma_sample] = draw_spiked_logitnormal(...
            [parameter_current(:).alpha_s],[parameter_current(:).beta_s],delay_sigma_bound);
%         t1e=toc;        time_rec(1)=time_rec(1)+t1e-t1s;
        
        shape_sample = randsample(1:length(z_dictionary),n_cell,true);
    
    
        delay_mu_sample_mat(:,s)=delay_mu_sample;
        delay_sigma_sample_mat(:,s)=delay_sigma_sample;
        gamma_sample_mat(:,s)=gamma_sample;
        gain_sample_mat(:,s)=gain_sample;
        shape_sample_mat(:,s)=shape_sample;
%         ts=toc;
%         % Calculate log pdf of prior and variational dist.
        logprior(:,s)=calculate_logpdf_spiked_logitnormal(prior_params,...
            logit_gamma,gamma_sample,logit_gain,gain_sample,...
            logit_delay_mu,delay_mu_sample,logit_delay_sigma,delay_sigma_sample,...
             delay_mu_bound,delay_sigma_bound,gamma_bound,gain_bound,spike_indicator);

        logvariational(:,s)=calculate_logpdf_spiked_logitnormal(parameter_current,...
            logit_gamma,gamma_sample,logit_gain,gain_sample,...
            logit_delay_mu,delay_mu_sample,logit_delay_sigma,delay_sigma_sample,...
             delay_mu_bound,delay_sigma_bound,gamma_bound,gain_bound);  
%         te=toc;        time_rec(2)=time_rec(2)+te-ts;
        
%         ts=toc;
        [dqdp_logit(:,s),dqdalpha(:,s),dqdbeta(:,s)] = get_dvariational_dist(gamma_sample,logit_gamma,...
            [parameter_current(:).alpha],[parameter_current(:).beta],gamma_constants,[parameter_current(:).p_logit],spike_indicator);
        [~,dqdalpha_gain(:,s),dqdbeta_gain(:,s)] = get_dvariational_dist(gain_sample,logit_gain,...
            [parameter_current(:).alpha_gain],[parameter_current(:).beta_gain],gain_constants);
        [~,dqdalpha_m(:,s),dqdbeta_m(:,s)] = get_dvariational_dist(delay_mu_sample,logit_delay_mu,...
            [parameter_current(:).alpha_m],[parameter_current(:).beta_m],delay_mu_constants);
        [~,dqdalpha_s(:,s),dqdbeta_s(:,s)] = get_dvariational_dist(delay_sigma_sample,logit_delay_sigma,...
            [parameter_current(:).alpha_s],[parameter_current(:).beta_s],delay_sigma_constants);
%         te=toc;        time_rec(3)=time_rec(3)+te-ts;
  
    end
    
        % Only thing that needs to change if we use the new way to get
        % intensity: 
%         ts=toc;
    [loglklh] = update_likelihood_dev_shape(power_size,mpp,...
        gamma_sample_mat,gain_sample_mat,delay_mu_sample_mat,delay_sigma_sample_mat,shape_sample_mat,...
            minimum_stim_threshold,background_rate,relevant_trials,lklh_func,...
        spike_curves,Tmax,z_dictionary);
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
           
