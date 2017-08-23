function [parameter_history] = fit_full_model_vi(...
    stim_size, mpp, background_rate, ...
    prob_trace,    stim_grid,...
    stim_scale,eff_stim_threshold,gain_bound,...
    variational_params,prior_params,C_threshold,stim_threshold,...
    S,epsilon,eta_logit,eta_beta,maxit)

% stim_threshold = 10;
% gain_bound.up=0.03;
% gain_bound.low=0.01;


% mpp=mpp_temp;
% stim_size=stim_size(:,find(cells_history{last_iter}));

n_cell=size(stim_size,2);n_trial=size(stim_size,1);
n_grid=size(prob_trace,2);

for i_cell = 1:n_cell
    variational_params(i_cell).pi = 0.01;
    variational_params(i_cell).p_logit = log(variational_params(i_cell).pi/(1-variational_params(i_cell).pi));
    variational_params(i_cell).log_alpha = 0;
    variational_params(i_cell).log_beta = 0;
    variational_params(i_cell).log_alpha_gain = 0;
    variational_params(i_cell).log_beta_gain = 0;
end
prior_params.pi0= 0.01*ones(n_cell,1);
prior_params.alpha0= ones(n_cell,1);
prior_params.beta0 = ones(n_cell,1);
prior_params.alpha0_gain= ones(n_cell,1);
prior_params.beta0_gain = ones(n_cell,1);

    

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

%parameter_history=struct([]);
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
v_log_alpha=[variational_params(:).log_alpha]';
v_log_beta=[variational_params(:).log_beta]';
v_log_alpha_gain=[variational_params(:).log_alpha_gain]';
v_log_beta_gain=[variational_params(:).log_beta_gain]';


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

% pickout the identifiable cells
while (changes > epsilon & iter<maxit)
    % Update the variational parameters using gradient descents
    
    %     t1=toc;
    v_alpha = exp(v_log_alpha);v_beta = exp(v_log_beta);
    v_alpha_gain = exp(v_log_alpha_gain);v_beta_gain = exp(v_log_beta_gain);
    
    v_pi_old = v_pi;v_p_logit_old=v_p_logit;
    v_log_alpha_old = v_log_alpha;v_log_beta_old = v_log_beta;
    v_log_alpha_gain_old = v_log_alpha_gain;v_log_beta_gain_old = v_log_beta_gain;
    
    
    parameter_history.pi(:,iter)=v_pi;
    parameter_history.p_logit(:,iter)=v_p_logit;
    parameter_history.alpha(:,iter)=v_alpha;
    parameter_history.beta(:,iter)=v_beta;
    parameter_history.alpha_gain(:,iter)=v_alpha_gain;
    parameter_history.beta_gain(:,iter)=v_beta_gain;
    %     v_pi=temp_rec_pi(:,iter);
    %     v_p_logit=temp_rec_p_logit(:,iter);
    %     v_alpha=temp_rec_alpha(:,iter);
    %     v_beta=temp_rec_beta(:,iter);
    
    iter=iter+1;
    alphaplusbeta = v_alpha+v_beta;
    psi_minus_alpha=psi(alphaplusbeta)-psi(v_alpha);
    psi_minus_beta=psi(alphaplusbeta)-psi(v_beta);
    
    alphaplusbeta_gain = v_alpha_gain+v_beta_gain;
    psi_minus_alpha_gain=psi(alphaplusbeta_gain)-psi(v_alpha_gain);
    psi_minus_beta_gain=psi(alphaplusbeta_gain)-psi(v_beta_gain);
    
%     t2=toc;time_record(1)=time_record(1)+t2-t1;
    for s= 1:S
        % Draw samples from the variational distribution given the current
        % parameters
        gamma_spike = rand(n_cell,1) > v_pi;
        gamma_slab = betarnd(v_alpha,v_beta,[n_cell 1])*(1-C_threshold) +C_threshold;
        gamma_sample = gamma_spike.*gamma_slab;
        gain_sample = betarnd(v_alpha_gain,v_beta_gain,[n_cell 1])*(gain_bound.up-gain_bound.low) +gain_bound.low;
        
        %bound gamma from 1 to avoid singularity
        gamma_sample_mat(:,s)=min(gamma_sample, 0.999);
        gain_sample_mat(:,s)=gain_sample;
        
        % Calculate the prior probability given the current sample of gamma
        % & gain 
        logprior(:,s)=log(max(0.001,prior_params.pi0)).*(gamma_sample==0)+...
            log(max(0.001, 1-prior_params.pi0)).*(gamma_sample>0)+ ...
            (gamma_sample>0).*log( min(1000,max(0.0001,...
            betapdf((gamma_sample),prior_params.alpha0,prior_params.beta0))))+...
            (gamma_sample>0).*log( min(1000,max(0.0001,...
         betapdf((gain_sample-gain_bound.low)/(gain_bound.up-gain_bound.low),prior_params.alpha0_gain,prior_params.beta0_gain)...
            /(gain_bound.up-gain_bound.low))));    
        
        % Calculate the probability of the variational distribution given the
        % current sample of gamma & gain 
        logvariational(:,s)=log(max(0.001,v_pi)).*(gamma_sample==0)+...
            log(max(0.001, 1-v_pi)).*(gamma_sample>0)+ ...
            (gamma_sample>0).*log( min(1000,max(0.0001,...
            betapdf((gamma_sample-C_threshold)/(1-C_threshold) ,v_alpha,v_beta)/(1-C_threshold))))...
            +(gamma_sample>0).*log( min(1000,max(0.0001,...
            betapdf((gain_sample-gain_bound.low)/(gain_bound.up-gain_bound.low),v_alpha_gain,v_beta_gain)...
            /(gain_bound.up-gain_bound.low))));
        
        
        dqdp_logit(:,s)= (gamma_sample==0)./(1+exp(v_p_logit))-...
            (gamma_sample>0).*exp(v_p_logit)./(1+exp(v_p_logit));
        dqdalpha(:,s)= (gamma_sample>0).*(log(max(0.001,(gamma_sample-C_threshold)/(1-C_threshold) )) +...
            psi_minus_alpha).*v_alpha;
        dqdalpha(gamma_sample==0,s)=0;
        dqdbeta(:,s)= (gamma_sample>0).*(log(max(0.001,1-(gamma_sample-C_threshold)/(1-C_threshold))) + ...
            psi_minus_beta).*v_beta;
        dqdbeta(gamma_sample==0,s)=0;
        dqdalpha_gain(:,s)= (gamma_sample>0).*...
            (log(max(0.001,(gain_sample-gain_bound.low)/(gain_bound.up-gain_bound.low))) +...
            psi_minus_alpha_gain).*v_alpha_gain;
        dqdalpha_gain(gamma_sample==0,s)=0;
        dqdbeta_gain(:,s)= (gamma_sample>0).*...
            (log(max(0.001,1-(gain_sample-gain_bound.low)/(gain_bound.up-gain_bound.low))) + ...
            psi_minus_beta_gain).*v_beta_gain;
        dqdbeta_gain(gamma_sample==0,s)=0;
        
    end
    fprintf('Iteration: %d; Sample drawn; ',iter);
%     t3=toc;time_record(2)=time_record(2)+t3-t2;
    
    for s=1:S
%         t3p=toc;
        gamma_sample=gamma_sample_mat(:,s);
        gain_sample=gain_sample_mat(:,s);
        loglklh_vec = zeros(n_trial,1);
        
        for  i_trial = 1:n_trial
            
            
            effective_stim= [stim_size(i_trial,:)'].*gain_sample;
            stimulated_cells = find(effective_stim>eff_stim_threshold);
            effective_stim=effective_stim(stimulated_cells );
            stim_index=max(1,round(effective_stim*stim_scale));
            
            prob_this_trial= (gamma_sample(stimulated_cells)*ones(1,n_grid)).*prob_trace(stim_index,:);
            [loss]=  lif_glm_firstspike_loglikelihood_for_VI(mpp(i_trial),...
                prob_this_trial,background_rate);
            loglklh_vec(i_trial)=loss;
        end
%         t4=toc;time_record(3)=time_record(3)+t4-t3p;
   
%           t5=toc;time_record(4)=time_record(4)+t5-t4;
        
        for i_cell = 1:n_cell
%             relevant_n_trial=length(relevant_trials);
%             fprintf('Cell: %d; num trials: %d;\n',i_cell,relevant_n_trial);
            loglklh(i_cell,s)=sum(loglklh_vec(relevant_trials{i_cell}));
        end
%          t6=toc;time_record(5)=time_record(5)+t6-t5;
   
    end
    
    fprintf(' Likelihood calculated;');
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
        );
        %v_pi = v_pi+eta*mean(dELBOdpi,2);
        if iter < 20
          v_p_logit(i_cell) = v_p_logit(i_cell)+eta_logit*mean(f_p_logit(i_cell,:)-a_constant*h_p_logit(i_cell,:));
            v_log_alpha(i_cell) = v_log_alpha(i_cell)+eta_beta*mean(f_alpha(i_cell,:)-a_constant*h_alpha(i_cell,:));
            v_log_beta(i_cell) = v_log_beta(i_cell)+eta_beta*mean(f_beta(i_cell,:)-a_constant*h_beta(i_cell,:));
            v_log_alpha(i_cell) = v_log_alpha_gain(i_cell)+...
                eta_beta*mean(f_alpha_gain(i_cell,:)-a_constant*h_alpha_gain(i_cell,:));
            v_log_beta(i_cell) = v_log_beta_gain(i_cell)+...
                eta_beta*mean(f_beta_gain(i_cell,:)-a_constant*h_beta_gain(i_cell,:));
        else
          v_p_logit(i_cell) = v_p_logit(i_cell)+(eta_beta/sqrt(iter*log(iter)))*mean(f_p_logit(i_cell,:)-a_constant*h_p_logit(i_cell,:));
            v_log_alpha(i_cell) = v_log_alpha(i_cell)+(eta_beta/sqrt(iter*log(iter)))*mean(f_alpha(i_cell,:)-a_constant*h_alpha(i_cell,:));
            v_log_beta(i_cell) = v_log_beta(i_cell)+(eta_beta/sqrt(iter*log(iter)))*mean(f_beta(i_cell,:)-a_constant*h_beta(i_cell,:));
           v_log_alpha(i_cell) = v_log_alpha(i_cell)+...
               (eta_beta/sqrt(iter*log(iter)))*mean(f_alpha_gain(i_cell,:)-a_constant*h_alpha_gain(i_cell,:));
            v_log_beta(i_cell) = v_log_beta(i_cell)+...
                (eta_beta/sqrt(iter*log(iter)))*mean(f_beta_gain(i_cell,:)-a_constant*h_beta_gain(i_cell,:));
        
        end
    end
    fprintf('Gradients obtained;');
    
%     t7=toc;time_record(6)=time_record(6)+t7-t6;
    v_pi = exp(v_p_logit)./(1+exp(v_p_logit));
    v_log_alpha=min(log(100),max(log(1e-2),v_log_alpha));
    v_log_beta=min(log(100),max(log(1e-2),v_log_beta));
    v_log_alpha_gain=min(log(100),max(log(1e-2),v_log_alpha_gain));
    v_log_beta_gain=min(log(100),max(log(1e-2),v_log_beta_gain));
    
    v_alpha = exp(v_log_alpha);
    v_beta = exp(v_log_beta);
    v_alpha_gain = exp(v_log_alpha_gain);
    v_beta_gain = exp(v_log_beta_gain);
    
    % Calculate the stopping criteriar
    %changes=sqrt(sum(mean(dELBOdpi,2).^2)+sum(mean(dELBOdalpha,2).^2)+sum(mean(dELBOdbeta,2).^2));
%     mean_gamma= (1-v_pi).*(C_threshold+ (1-C_threshold)*v_alpha./(v_alpha+v_beta));
%     mean_gamma_old= (1-v_pi_old).*(C_threshold+ (1-C_threshold)*v_alpha_old./(v_alpha_old+v_beta_old));
    
    %changes = sqrt(sum((mean_gamma-mean_gamma_old).^2)/sum(mean_gamma_old.^2 ));
    changes=  sum(abs(v_pi_old-v_pi))+...
        sum(abs(v_log_alpha_old-v_log_alpha))+...
        sum(abs(v_log_beta_old-v_log_beta)) +...
    sum(abs(v_log_alpha_gain_old-v_log_alpha))+...
        sum(abs(v_log_beta_gain_old-v_log_beta_gain));
    changes= changes/(sum(abs(v_log_alpha_old))+sum(abs(v_log_beta_old))+sum(abs(v_pi_old))+...
        sum(abs(v_log_alpha_gain_old))+sum(abs(v_log_beta_gain_old))); 
    
    change_history(iter)=changes;
    
    fprintf('Change: %d;\n',changes)
end