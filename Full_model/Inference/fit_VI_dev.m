% function [parameter_history] = fit_VI_dev(...
%     stim_size, mpp, background_rate, ...
%     variational_params,prior_params,...
%     inference_params,prior_info,varargin)
%%
% stim_size=designs_remained;
% background_rate=experiment_setup.patched_neuron.background_rate;
mpp=mpp_remained;
inference_params=newmodel_profile.inference_params;
%%

% if ~isempty(varargin) && ~isempty(varargin{1})
%    spike_indicator= varargin{1};
% else
%     spike_indicator= false;
% end

% No need to match these scales anymore:
% intensity_grid=prior_info.induced_intensity.intensity_grid;
% stim_scale=prior_info.induced_intensity.stim_scale;
% minimum_stim_threshold=prior_info.induced_intensity.minimum_stim_threshold;

epsilon=inference_params.convergence_threshold;
S=inference_params.MCsamples_for_gradient;
eta=inference_params.step_size;
eta_max=inference_params.step_size_max;
maxit=inference_params.maxit;
lklh_func=inference_params.likelihood;
%% Undefined parameters: 
n_grid=300;


dqdp_logit=zeros(n_cell,S);
dqdalpha=zeros(n_cell,S);dqdbeta=zeros(n_cell,S);
dqdalpha_gain=zeros(n_cell,S);dqdbeta_gain=zeros(n_cell,S);
%%
%stim_size=designs_remained; mpp= mpp_remained;
n_cell=size(stim_size,2);n_trial=size(stim_size,1);
% initialize storages 
sum_of_logs=zeros(S,1);logvariational=zeros(n_cell,S);
logprior=zeros(n_cell,S);loglklh=zeros(n_cell,S);
change_history=epsilon+1;

param_names=fieldnames(variational_params(1));
parameter_sample=struct;
for i_param = 1:length(param_names)
    parameter_sample.(param_names{i_param})=zeros(n_cell,S);
end

%
parameter_history=variational_params;
parameter_current=variational_params;

% % find relevant trials for each cell 
% relevant_trials = cell(n_cell,1);
% for i_cell = 1:n_cell
%     relevant_trials{i_cell}=find(stim_size(:,i_cell)>(minimum_stim_threshold/gain_bound.up));
% end
changes=1;iteration = 1;
%%
% while (change_history(iteration) > epsilon && iteration<maxit)
    parameter_history(iteration,:)=parameter_current;
    iteration=iteration+1;
    
    % precalculate some quantities 
    [gamma_constants] = get_logitnormal_constants([parameter_current(:).alpha],[parameter_current(:).beta]);
    [gain_constants] = get_logitnormal_constants([parameter_current(:).alpha_gain],[parameter_current(:).beta_gain]);
    %
    
    % Output the transformed parameters for evaluating the likelihood
    % Also output the untransformed parameters for evaluating variational
    % probability 
    [parameter_sample, raw_sample]=draw_variational_samples(parameter_sample,parameter_current,inference_params.bounds);
    
    % Evaluating the variational densities and their derivatives 
    % Revise this part so that it is only one function
    
    get_prior_distribution();
    
      get_variational_distribution(parameter_current,...
            parameter_sample,raw_sample);
   %this time with derivatives 
    
    for s= 1:S
        % drawing samples 
        % Calculate log pdf of prior and variational dist.
        logprior(:,s)=calculate_logpdf_spiked_logitnormal(prior_params,...
            logit_gamma,gamma_sample,logit_gain,gain_sample,...
            gamma_bound,gain_bound,spike_indicator);
        logvariational(:,s)=calculate_logpdf_spiked_logitnormal(parameter_current,...
            logit_gamma,gamma_sample,logit_gain,gain_sample,...
            gamma_bound,gain_bound);
        [dqdp_logit(:,s),dqdalpha(:,s),dqdbeta(:,s)] = get_dvariational_dist(gamma_sample,logit_gamma,...
            [parameter_current(:).alpha],[parameter_current(:).beta],gamma_constants,[parameter_current(:).p_logit],spike_indicator);
        [~,dqdalpha_gain(:,s),dqdbeta_gain(:,s)] = get_dvariational_dist(gain_sample,logit_gain,...
            [parameter_current(:).alpha_gain],[parameter_current(:).beta_gain],gain_constants);
    end
    
    %%
%     [mean(gamma_sample_mat)  [parameter_current(:).alpha]]
%     [var(gamma_sample_mat)  exp(2*[parameter_current(:).beta])]
%     [mean(gain_sample_mat)  [parameter_current(:).alpha_gain]]
%     [var(gain_sample_mat)  exp(2*[parameter_current(:).beta_gain])]
    
    [loglklh] = update_likelihood(gamma_sample_mat,gain_sample_mat,stim_size,mpp,...
        intensity_grid,minimum_stim_threshold,stim_scale,background_rate,relevant_trials,lklh_func);

    [parameter_current, change_history(iteration)]= update_parameters_logitnormal(...
        parameter_current,loglklh, logprior, logvariational,...
        dqdp_logit, dqdalpha, dqdbeta, dqdalpha_gain, dqdbeta_gain,...
        iteration,eta,eta_max,spike_indicator,spike_indicator);
% parameter_current

% fprintf('Iteration %d; change %d;\n',iteration, change_history(iteration))
% end
% fprintf('VI fitted after %d iteration;\n',iteration)
           
