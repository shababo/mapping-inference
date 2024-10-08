function [parameter_history] = fit_working_model_vi_spike(...
    stim_size, mpp, background_rate, ...
    prob_trace_full, stim_grid,...
    stim_scale,eff_stim_threshold,gain_bound,...
    variational_params,prior_params,C_threshold,stim_threshold,...
    stim_size_neighbours,gamma_neighbours,gain_neighbours,...
    S,epsilon,eta_beta,eta_max,maxit,lklh_func)

% stim_size=designs_remained; mpp=mpp_remained;

%    stim_size_neighbours= designs_neighbours;
%  stim_size;
%  mpp=mpp_temp;

% mpp=mpp_temp;
% stim_size=stim_size(:,find(cells_history{last_iter}));
%%
n_cell=size(stim_size,2);n_trial=size(stim_size,1);
n_grid=size(prob_trace_full,2);

sum_of_logs=zeros(S,1);logvariational=zeros(n_cell,S);
logprior=zeros(n_cell,S);loglklh=zeros(n_cell,S);

dqdp_logit=zeros(n_cell,S);
dqdalpha=zeros(n_cell,S);dqdbeta=zeros(n_cell,S);
dqdalpha_gain=zeros(n_cell,S);dqdbeta_gain=zeros(n_cell,S);

f_p_logit=zeros(n_cell,S);
f_alpha=zeros(n_cell,S);f_beta=zeros(n_cell,S);
f_alpha_gain=zeros(n_cell,S);f_beta_gain=zeros(n_cell,S);
h_p_logit=zeros(n_cell,S);
h_alpha=zeros(n_cell,S);h_beta=zeros(n_cell,S);
h_alpha_gain=zeros(n_cell,S);h_beta_gain=zeros(n_cell,S);

% parameter_history=struct([]);
parameter_history.pi=zeros(n_cell,1);
parameter_history.p_logit=zeros(n_cell,1);
parameter_history.alpha=zeros(n_cell,1);
parameter_history.beta=zeros(n_cell,1);
parameter_history.alpha_gain=zeros(n_cell,1);
parameter_history.beta_gain=zeros(n_cell,1);

change_history=zeros(1,1);

gamma_sample_mat=zeros(n_cell,S);
gain_sample_mat=zeros(n_cell,S);



v_pi =[variational_params(:).pi]';
v_p_logit =[variational_params(:).p_logit]';
v_log_alpha=[variational_params(:).alpha]';
v_log_beta=log([variational_params(:).beta]');
v_log_alpha_gain=[variational_params(:).alpha_gain]';
v_log_beta_gain=log([variational_params(:).beta_gain]');


% tic;
%
% time_record = zeros(10,1);

%------------------------------%
% Store the relevant trials for each cell
relevant_trials = cell(n_cell,1);
for i_cell = 1:n_cell
    relevant_trials{i_cell}=find(stim_size(:,i_cell)>stim_threshold);
end
changes=1;iter = 1;
%
% pickout the identifiable cells
while (changes > epsilon & iter<maxit)
    % Update the variational parameters using gradient descents
    v_alpha = v_log_alpha;v_beta = exp(v_log_beta);
    v_alpha_gain = v_log_alpha_gain;v_beta_gain = exp(v_log_beta_gain);
      v_pi_old = v_pi;v_p_logit_old=v_p_logit;
    v_log_alpha_old = v_log_alpha;v_log_beta_old = v_log_beta;
    v_log_alpha_gain_old = v_log_alpha_gain;v_log_beta_gain_old = v_log_beta_gain;
    
        parameter_history.pi(:,iter)=v_pi;
    parameter_history.p_logit(:,iter)=v_p_logit;
    parameter_history.alpha(:,iter)=v_alpha;
    parameter_history.beta(:,iter)=v_beta;
    parameter_history.alpha_gain(:,iter)=v_alpha_gain;
    parameter_history.beta_gain(:,iter)=v_beta_gain;
    
    iter=iter+1;
    mu_by_sigma2 = v_alpha./(v_beta.^2);
    sigma_inv= 1./(v_beta);
    sigma_inv2= 1./(v_beta.^2);
    sigma_inv3= 1./(v_beta.^3);
    
    mu_by_sigma2_gain = v_alpha_gain./(v_beta_gain.^2);
    sigma_inv_gain= 1./(v_beta_gain);
    sigma_inv2_gain= 1./(v_beta_gain.^2);
    sigma_inv3_gain= 1./(v_beta_gain.^3);
    
    for s= 1:S
        % Draw samples from the variational distribution given the current
        % parameters
       logit_gamma=normrnd(v_alpha,v_beta,[n_cell 1]);
          gamma_spike = rand(n_cell,1) > v_pi;
        gamma_slab = exp(logit_gamma)./(1+exp(logit_gamma))*(1-C_threshold) +C_threshold;
        gamma_sample = gamma_slab.*gamma_spike;
       logit_gain=normrnd(v_alpha_gain,v_beta_gain,[n_cell 1]);
        gain_sample = exp( logit_gain)./(1+exp( logit_gain))*(gain_bound.up-gain_bound.low) +gain_bound.low;
        % to avoid singularity
        gain_sample = max(gain_sample,gain_bound.low+0.001);
        gain_sample = min(gain_sample,gain_bound.up-0.001);
        %bound gamma from 1 to avoid singularity
        gamma_sample_mat(:,s)=min(gamma_sample, 0.999);
        gain_sample_mat(:,s)=gain_sample;
        % Calculate the prior probability given the current sample of gamma
        % & gain
        
        % Calculate the probability of the variational distribution given the
        % current sample of gamma & gain
        logprior(:,s)=log(max(0.001,prior_params.pi0)).*(gamma_sample==0)+...
            log(max(0.001, 1-prior_params.pi0)).*(gamma_sample>0)+ ...
            (gamma_sample>0).*log( min(1000,max(0.0001,...
            normpdf(logit_gamma,prior_params.alpha0,prior_params.beta0)./(gamma_sample.*(1-gamma_sample))/(1-C_threshold)  )))...
            +(gamma_sample>0).*log( min(1000,max(0.0001,...
            normpdf(logit_gain,prior_params.alpha0_gain,prior_params.beta0_gain)...
            ./( (gain_sample-gain_bound.low ).*(gain_bound.up-gain_sample)) *(gain_bound.up-gain_bound.low))));
        %-----------------------------------------%
        % Need to change the following functions when changing the
        % variational families :
          
        logvariational(:,s)=log(max(0.001,v_pi)).*(gamma_sample==0)+...
        log(max(0.001, 1-v_pi)).*(gamma_sample>0).*log( min(1000,max(0.0001,...
            normpdf(logit_gamma,v_alpha,v_beta)./(gamma_sample.*(1-gamma_sample))/(1-C_threshold)  )))...
            +log( min(1000,max(0.0001,...
            normpdf(logit_gain,v_alpha_gain,v_beta_gain)...
            ./( (gain_sample-gain_bound.low ).*(gain_bound.up-gain_sample)) *(gain_bound.up-gain_bound.low))));
        
                dqdp_logit(:,s)= (gamma_sample==0)./(1+exp(v_p_logit))-...
            (gamma_sample>0).*exp(v_p_logit)./(1+exp(v_p_logit));

        dqdalpha(:,s)= (-mu_by_sigma2+logit_gamma.*sigma_inv2);
                dqdalpha(gamma_sample==0,s)=0;
        dqdbeta(:,s)= (-sigma_inv+sigma_inv3.*(logit_gamma-v_alpha).^2).*v_beta;
                dqdbeta(gamma_sample==0,s)=0;
        
        dqdalpha_gain(:,s)= (-mu_by_sigma2_gain+logit_gain.*sigma_inv2_gain);
                dqdalpha_gain(gamma_sample==0,s)=0;
        dqdbeta_gain(:,s)= (-sigma_inv_gain+sigma_inv3_gain.*(logit_gain-v_alpha_gain).^2).*v_beta_gain;
                dqdbeta_gain(gamma_sample==0,s)=0;
        %--------------------------------------------------------%
        
    end
    %     fprintf('Iteration: %d; Sample drawn; ',iter);
    %     t3=toc;time_record(2)=time_record(2)+t3-t2;
    
    for s=1:S
        %         t3p=toc;
        gamma_sample=gamma_sample_mat(:,s);
        gain_sample=gain_sample_mat(:,s);
        
        % Need to account for the neighbour cells that are not in this
        % group:
        loglklh_vec = zeros(n_trial,1);
        for  i_trial = 1:n_trial
            if ~isempty(gamma_neighbours)
                stim_temp =[stim_size(i_trial,:)';stim_size_neighbours(i_trial,:)'];
                gain_temp=[gain_sample;gain_neighbours];
                gamma_temp=[gamma_sample;gamma_neighbours];
            else
                stim_temp =stim_size(i_trial,:)';
                gain_temp=gain_sample;
                gamma_temp=gamma_sample;
            end
            effective_stim= stim_temp.*gain_temp;
            stimulated_cells = find(effective_stim>eff_stim_threshold/2);
            effective_stim=effective_stim(stimulated_cells );
            stim_index=min(size(prob_trace_full,1), max(1,round(effective_stim*stim_scale)));
            if ~isempty(stimulated_cells)
                prob_this_trial= prob_trace_full(stim_index,:);
                prob_this_trial=[background_rate*ones(1, size(prob_this_trial,2)); prob_this_trial];
            else
                prob_this_trial=[background_rate*ones(1, size(prob_trace_full,2))];
            end
            prob_collapsed=sum(prob_this_trial,2);
            n_events = length(mpp(i_trial).times);
            [lklh]= lklh_func(n_events,...
                [1;gamma_temp(stimulated_cells)],prob_collapsed);
            
            loglklh_vec(i_trial)=log(lklh);
        end
        for i_cell = 1:n_cell
            loglklh(i_cell,s)=sum(loglklh_vec(relevant_trials{i_cell}));
        end
        %          t6=toc;time_record(5)=time_record(5)+t6-t5;
        
    end
    
    %     fprintf(' Likelihood calculated;');
    for i_cell = 1:n_cell
        sum_of_logs = loglklh(i_cell,:)+logprior(i_cell,:)-logvariational(i_cell,:);
        
          f_p_logit(i_cell,:)=sum_of_logs.*dqdp_logit(i_cell,:);
        f_alpha(i_cell,:)= sum_of_logs.*dqdalpha(i_cell,:);
        f_beta(i_cell,:)= sum_of_logs.*dqdbeta(i_cell,:);
        f_alpha_gain(i_cell,:)= sum_of_logs.*dqdalpha_gain(i_cell,:);
        f_beta_gain(i_cell,:)= sum_of_logs.*dqdbeta_gain(i_cell,:);
        
           h_p_logit(i_cell,:)=dqdp_logit(i_cell,:);
        h_alpha(i_cell,:)= dqdalpha(i_cell,:);
        h_beta(i_cell,:)= dqdbeta(i_cell,:);
        h_alpha_gain(i_cell,:)= dqdalpha_gain(i_cell,:);
        h_beta_gain(i_cell,:)= dqdbeta_gain(i_cell,:);
        
        % Calculate the constant a
        a_constant = (quick_cov(f_p_logit(i_cell,:),h_p_logit(i_cell,:))+...
            quick_cov(f_alpha(i_cell,:),h_alpha(i_cell,:))+...
            quick_cov(f_beta(i_cell,:),h_beta(i_cell,:))+...
            quick_cov(f_alpha_gain(i_cell,:),h_alpha_gain(i_cell,:))+...
            quick_cov(f_beta_gain(i_cell,:),h_beta_gain(i_cell,:))...
            )/...
            (quick_cov(h_p_logit(i_cell,:),h_p_logit(i_cell,:))+...
            quick_cov(h_alpha(i_cell,:),h_alpha(i_cell,:))+...
            quick_cov(h_beta(i_cell,:),h_beta(i_cell,:))+...
            quick_cov(h_alpha_gain(i_cell,:),h_alpha_gain(i_cell,:))+...
            quick_cov(h_beta_gain(i_cell,:),h_beta_gain(i_cell,:))...
            +0.01); % to prevent singularity 
        %v_pi = v_pi+eta*mean(dELBOdpi,2);
        grad_logit=(eta_beta/sqrt(iter*log(iter)))*mean(f_p_logit(i_cell,:)-a_constant*h_p_logit(i_cell,:));
        grad_log_alpha = (eta_beta/sqrt(iter*log(iter)))*mean(f_alpha(i_cell,:)-a_constant*h_alpha(i_cell,:));
       grad_log_beta = (eta_beta/sqrt(iter*log(iter)))*mean(f_beta(i_cell,:)-a_constant*h_beta(i_cell,:));
       grad_log_alpha_gain=(eta_beta/sqrt(iter*log(iter)))*mean(f_alpha_gain(i_cell,:)-a_constant*h_alpha_gain(i_cell,:));
       grad_log_beta_gain=(eta_beta/sqrt(iter*log(iter)))*mean(f_beta_gain(i_cell,:)-a_constant*h_beta_gain(i_cell,:));
       grad_max = max(abs([grad_logit grad_log_alpha grad_log_beta grad_log_alpha_gain grad_log_beta_gain]));
   
        if grad_max > eta_max
           grad_scale= grad_max/eta_max;
            grad_logit= grad_logit/grad_scale;
           grad_log_alpha = grad_log_alpha/grad_scale;
           grad_log_beta = grad_log_beta/grad_scale;
           grad_log_alpha_gain=grad_log_alpha_gain/grad_scale;
           grad_log_beta_gain=grad_log_beta_gain/grad_scale;
       end
              v_p_logit(i_cell) = v_p_logit(i_cell)+grad_logit;
             v_log_alpha(i_cell) = v_log_alpha(i_cell)+grad_log_alpha;
            v_log_beta(i_cell) = v_log_beta(i_cell)+grad_log_beta;
           v_log_alpha_gain(i_cell) = v_log_alpha_gain(i_cell)+...
               grad_log_alpha_gain;
            v_log_beta_gain(i_cell) = v_log_beta_gain(i_cell)+...
                grad_log_beta_gain;
    end
    %     fprintf('Gradients obtained;');
    
       v_pi = exp(v_p_logit)./(1+exp(v_p_logit));
    v_log_beta=min(log(100),max(log(1e-2),v_log_beta));
    v_log_beta_gain=min(log(100),max(log(1e-2),v_log_beta_gain));
    v_beta = exp(v_log_beta);
    v_beta_gain = exp(v_log_beta_gain);
    v_alpha = v_log_alpha;
    v_alpha_gain = v_log_alpha_gain;
    
    % Calculate the stopping criteriar
    %changes=sqrt(sum(mean(dELBOdpi,2).^2)+sum(mean(dELBOdalpha,2).^2)+sum(mean(dELBOdbeta,2).^2));
    %  mean_gamma= (1-v_pi).*(C_threshold+ (1-C_threshold)*v_alpha./(v_alpha+v_beta));
    %  mean_gamma=[1 1]';
    %     mean_gamma_old= (1-v_pi_old).*(C_threshold+ (1-C_threshold)*v_alpha_old./(v_alpha_old+v_beta_old));
    
    %changes = sqrt(sum((mean_gamma-mean_gamma_old).^2)/sum(mean_gamma_old.^2 ));
    changes=  sum(abs(v_log_alpha_old-v_log_alpha))+...
        sum(abs(v_log_beta_old-v_log_beta)) +...
        sum(abs(v_log_alpha_gain_old-v_log_alpha_gain))+...
        sum(abs(v_log_beta_gain_old-v_log_beta_gain));
    changes= changes/(sum(abs(v_log_alpha_old))+sum(abs(v_log_beta_old))+...
        sum(abs(v_log_alpha_gain_old))+sum(abs(v_log_beta_gain_old)));
    
    change_history(iter)=changes;
    
        %fprintf('Change: %d;\n',changes)
end
fprintf('VI fitted after %d iteration;\n',iter)
end

%% Debug section:
%  [mean_gamma_temp, var_gamma_temp] = calculate_posterior_mean(...
%         parameter_history.alpha(:,end),parameter_history.beta(:,end),0,1);
%     [mean_gain_temp, var_gain_temp] = calculate_posterior_mean(...
%         parameter_history.alpha_gain(:,end),parameter_history.beta_gain(:,end),gain_bound.low,gain_bound.up);
%    
% %
% gamma_truth
% mean_gamma_temp
% gain_truth
% mean_gain_temp
