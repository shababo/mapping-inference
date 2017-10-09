function [parameter_history,change_history] = fit_working_model_vi(...
    designs,outputs, background_rate, ...
    variational_params,prior_params,C_threshold,...
    designs_neighbours,gamma_neighbours,...
    S,epsilon,eta_logit,eta_beta,maxit,lklh_func)

vf_type =2; % use logit-normal distribution 

%   designs=designs_remained;
%   outputs=outputs_remained;
%   background_rate=background_rt;
% lklh_func=@calculate_likelihood_sum_bernoulli;
n_cell=size(designs,2);
n_trial=size(designs,1);


sum_of_logs=zeros(S,1);
logvariational=zeros(n_cell,S);
logprior=zeros(n_cell,S);
loglklh=zeros(n_cell,S);

dqdp_logit=zeros(n_cell,S);
dqdalpha=zeros(n_cell,S);
dqdbeta=zeros(n_cell,S);

f_p_logit=zeros(n_cell,S);
f_alpha=zeros(n_cell,S);
f_beta=zeros(n_cell,S);

h_p_logit=zeros(n_cell,S);
h_dalpha=zeros(n_cell,S);
h_beta=zeros(n_cell,S);

%parameter_history=struct([]);
% parameter_history.pi=zeros(n_cell,1);
% parameter_history.p_logit=zeros(n_cell,1);
parameter_history.alpha=zeros(n_cell,1);
parameter_history.beta=zeros(n_cell,1);
change_history=zeros(1,1);

gamma_sample_mat=zeros(n_cell,S);

if vf_type ==1
%     v_pi =[variational_params(:).pi]';
%     v_p_logit =[variational_params(:).p_logit]';
    v_log_alpha=log([variational_params(:).alpha]');
    v_log_beta=log([variational_params(:).beta]');
elseif vf_type == 2
%     v_pi =[variational_params(:).pi]';
%     v_p_logit =[variational_params(:).p_logit]';
    v_log_alpha=[variational_params(:).alpha]';
    v_log_beta=log([variational_params(:).beta]');
end

if isempty(designs_neighbours)
    data_matrix = [outputs designs];
else
    data_matrix = [outputs designs designs_neighbours];
end

[data_unique,~,data_index] = unique(data_matrix, 'rows');
n_unique = size(data_unique,1);
loglikelihood_unique = zeros(n_unique,1);


relevant_trials = cell(n_cell,1);
for i_cell = 1:n_cell
    relevant_trials{i_cell}=find(designs(:,i_cell)>0);
end
% tic;
%
% time_record = zeros(10,1);
changes=1;

iter = 1;
% pickout the identifiable cells
while (changes > epsilon & iter<maxit)
    % Update the variational parameters using gradient descents
    
    %     t1=toc;
    
    if vf_type == 1
        v_alpha = exp(v_log_alpha);
        v_beta = exp(v_log_beta);
%         v_pi_old = v_pi;v_p_logit_old=v_p_logit;
        v_log_alpha_old = v_log_alpha;v_log_beta_old = v_log_beta;
        v_alpha_old=v_alpha;v_beta_old=v_beta;
    elseif vf_type == 2
        v_alpha = v_log_alpha;
        v_beta = exp(v_log_beta);
%         v_pi_old = v_pi;v_p_logit_old=v_p_logit;
        v_log_alpha_old = v_log_alpha;v_log_beta_old = v_log_beta;
        v_alpha_old=v_alpha;v_beta_old=v_beta; 
    end
    
%     parameter_history.pi(:,iter)=v_pi;
%     parameter_history.p_logit(:,iter)=v_p_logit;
    parameter_history.alpha(:,iter)=v_alpha;
    parameter_history.beta(:,iter)=v_beta;
    
    %     v_pi=temp_rec_pi(:,iter);
    %     v_p_logit=temp_rec_p_logit(:,iter);
    %     v_alpha=temp_rec_alpha(:,iter);
    %     v_beta=temp_rec_beta(:,iter);
    iter=iter+1;
    
     if vf_type== 1 % beta 
    alphaplusbeta = v_alpha+v_beta;
    psi_minus_alpha=psi(alphaplusbeta)-psi(v_alpha);
    psi_minus_beta=psi(alphaplusbeta)-psi(v_beta);
    elseif vf_type==2 % logit-normal
        mu_by_sigma2 = v_alpha./(v_beta.^2);
        sigma_inv= 1./(v_beta);
        sigma_inv2= 1./(v_beta.^2);
        sigma_inv3= 1./(v_beta.^3);
     end
    
    %     t2=toc;time_record(1)=time_record(1)+t2-t1;
    for s= 1:S
        % Draw samples from the variational distribution given the current
        % parameters
        if vf_type == 1
%             gamma_spike = rand(n_cell,1) > v_pi;
            gamma_slab = betarnd(v_alpha,v_beta,[n_cell 1])*(1-C_threshold) +C_threshold;
            gamma_sample = gamma_slab;
        elseif vf_type == 2
%             gamma_spike = rand(n_cell,1) > v_pi;
            temp=normrnd(v_alpha,v_beta,[n_cell 1]);
            gamma_slab = exp(temp)./(1+exp(temp))*(1-C_threshold) +C_threshold;
            gamma_sample = gamma_slab;
        end
        
        %bound gamma from 1 to avoid singularity
        gamma_sample_mat(:,s)=min(gamma_sample, 0.999);
        % Calculate the prior probability given the current sample of gamma
        if vf_type==1
            logprior(:,s)=log(max(0.001,prior_params.pi0)).*(gamma_sample==0)+...
                log(max(0.001, 1-prior_params.pi0)).*(gamma_sample>0)+ ...
                (gamma_sample>0).*log( min(1000,max(0.0001,...
                betapdf((gamma_sample),prior_params.alpha0,prior_params.beta0 ))));
        elseif vf_type==2
            logit_gamma = log( 1./ ( (1-C_threshold)./(gamma_sample-C_threshold)-1 ));
            logprior(:,s)=log(max(0.001,prior_params.pi0)).*(gamma_sample==0)+...
                log(max(0.001, 1-prior_params.pi0)).*(gamma_sample>0)+ ...
                (gamma_sample>0).*log( min(1000,max(0.0001,...
                normpdf(logit_gamma,prior_params.alpha0,prior_params.beta0)./(gamma_sample.*(1-gamma_sample))/(1-C_threshold)  )));
            
        end
        %          logprior(:,s)=log(prior_params.pi0).*(gamma_sample==0)+log(1-prior_params.pi0).*(gamma_sample>0)+ ...
        %             (gamma_sample>0).*log(betapdf(gamma_sample,prior_params.alpha0,prior_params.beta0));
        % Calculate the probability of the variational distribution given the
        % current sample of gamma
        if vf_type == 1
        logvariational(:,s)=  (gamma_sample>0).*log( min(1000,max(0.0001,...
            betapdf((gamma_sample-C_threshold)/(1-C_threshold) ,v_alpha,v_beta)/(1-C_threshold))));
        
%             dqdp_logit(:,s)= (gamma_sample==0)./(1+exp(v_p_logit))-...
%             (gamma_sample>0).*exp(v_p_logit)./(1+exp(v_p_logit));
        dqdalpha(:,s)= (gamma_sample>0).*(log(max(0.001,(gamma_sample-C_threshold)/(1-C_threshold) )) +...
            psi_minus_alpha).*v_alpha;
        dqdalpha(gamma_sample==0,s)=0;
        dqdbeta(:,s)= (gamma_sample>0).*(log(max(0.001,1-(gamma_sample-C_threshold)/(1-C_threshold))) + ...
            psi_minus_beta).*v_beta;
        dqdbeta(gamma_sample==0,s)=0;
        
        elseif vf_type == 2 % logit-normal distribution 
            
        % Calculate the probability of the variational distribution given the
        % current sample of gamma & gain 
%         logit_gamma = log( 1./ ( (1-C_threshold)./(gamma_sample-C_threshold)-1 ));
        logvariational(:,s)=   (gamma_sample>0).*log( min(1000,max(0.0001,...
           normpdf(logit_gamma,v_alpha,v_beta)./(gamma_sample.*(1-gamma_sample))/(1-C_threshold)  )));
        
%         dqdp_logit(:,s)= (gamma_sample==0)./(1+exp(v_p_logit))-...
%             (gamma_sample>0).*exp(v_p_logit)./(1+exp(v_p_logit));
        
        dqdalpha(:,s)= (gamma_sample>0).*(-mu_by_sigma2+logit_gamma.*sigma_inv2);
        dqdalpha(gamma_sample==0,s)=0;
        dqdbeta(:,s)= (gamma_sample>0).*(-sigma_inv+sigma_inv3.*(logit_gamma-v_alpha).^2).*v_beta;
        dqdbeta(gamma_sample==0,s)=0;
        end
    end
    fprintf('Iteration: %d; Sample drawn; ',iter);
    %     t3=toc;time_record(2)=time_record(2)+t3-t2;
    
    for s=1:S
        %         t3p=toc;
        if isempty(gamma_neighbours)
            gamma_temp=gamma_sample_mat(:,s);
        else
            gamma_temp=[gamma_sample_mat(:,s); gamma_neigbours];
        end
        for i_data = 1:n_unique
            n_events=data_unique(i_data,1);
            [lklh]= lklh_func(n_events,...
                [1;gamma_temp],[background_rate data_unique(i_data,2:end)]');
            loglikelihood_unique(i_data)=log(lklh);
        end
        %         t4=toc;time_record(3)=time_record(3)+t4-t3p;
        
        % Calculate the loglikelihood given the current samples of gamma
        % Need to change how we estimate the likelihood
        loglklh_vec = zeros(n_trial,1);
        for i_trial = 1:n_trial
            loglklh_vec(i_trial)=loglikelihood_unique(data_index(i_trial));
        end
        %         t5=toc;time_record(4)=time_record(4)+t5-t4;
        
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
        
%         f_p_logit(i_cell,:)=sum_of_logs.*dqdp_logit(i_cell,:);
        f_alpha(i_cell,:)= sum_of_logs.*dqdalpha(i_cell,:);
        f_beta(i_cell,:)= sum_of_logs.*dqdbeta(i_cell,:);
        
%         h_p_logit(i_cell,:)=dqdp_logit(i_cell,:);
        h_alpha(i_cell,:)= dqdalpha(i_cell,:);
        h_beta(i_cell,:)= dqdbeta(i_cell,:);
        
        % Calculate the constant a
        a_constant = (quick_cov(f_alpha(i_cell,:),h_alpha(i_cell,:))+quick_cov(f_beta(i_cell,:),h_beta(i_cell,:)))/...
            (quick_cov(h_alpha(i_cell,:),h_alpha(i_cell,:))+...
            quick_cov(h_beta(i_cell,:),h_beta(i_cell,:)));
        %v_pi = v_pi+eta*mean(dELBOdpi,2);
        if iter < 20
%             v_p_logit(i_cell) = v_p_logit(i_cell)+eta_logit*mean(f_p_logit(i_cell,:)-a_constant*h_p_logit(i_cell,:));
            v_log_alpha(i_cell) = v_log_alpha(i_cell)+eta_beta*mean(f_alpha(i_cell,:)-a_constant*h_alpha(i_cell,:));
            v_log_beta(i_cell) = v_log_beta(i_cell)+eta_beta*mean(f_beta(i_cell,:)-a_constant*h_beta(i_cell,:));
        else
%             v_p_logit(i_cell) = v_p_logit(i_cell)+(eta_beta/sqrt(iter*log(iter)))*mean(f_p_logit(i_cell,:)-a_constant*h_p_logit(i_cell,:));
            v_log_alpha(i_cell) = v_log_alpha(i_cell)+(eta_beta/sqrt(iter*log(iter)))*mean(f_alpha(i_cell,:)-a_constant*h_alpha(i_cell,:));
            v_log_beta(i_cell) = v_log_beta(i_cell)+(eta_beta/sqrt(iter*log(iter)))*mean(f_beta(i_cell,:)-a_constant*h_beta(i_cell,:));
        end
    end
    fprintf('Gradients obtained;');
    
    %     t7=toc;time_record(6)=time_record(6)+t7-t6;
    
    if vf_type == 1
%     v_pi = exp(v_p_logit)./(1+exp(v_p_logit));
    v_log_alpha=min(log(100),max(log(1e-2),v_log_alpha));
    v_log_beta=min(log(100),max(log(1e-2),v_log_beta));
    
    v_alpha = exp(v_log_alpha);
    v_beta = exp(v_log_beta);
   elseif vf_type == 2
    
%     v_pi = exp(v_p_logit)./(1+exp(v_p_logit));
    v_log_beta=min(log(100),max(log(1e-2),v_log_beta));
    
    v_beta = exp(v_log_beta);
    v_alpha = v_log_alpha;
   
    end

    
    % Calculate the stopping criteriar
    %changes=sqrt(sum(mean(dELBOdpi,2).^2)+sum(mean(dELBOdalpha,2).^2)+sum(mean(dELBOdbeta,2).^2));
    %mean_gamma= (1-v_pi).*(C_threshold+ (1-C_threshold)*v_alpha./(v_alpha+v_beta));
    %     mean_gamma_old= (1-v_pi_old).*(C_threshold+ (1-C_threshold)*v_alpha_old./(v_alpha_old+v_beta_old));
    
    %changes = sqrt(sum((mean_gamma-mean_gamma_old).^2)/sum(mean_gamma_old.^2 ));
    changes=          sum(abs(v_log_alpha_old-v_log_alpha))+...
        sum(abs(v_log_beta_old-v_log_beta)) ;
    changes= changes/(sum(abs(v_log_alpha_old))+sum(abs(v_log_beta_old)));
    
    change_history(iter)=changes;
    
    fprintf('Change: %d;\n',changes)
end
