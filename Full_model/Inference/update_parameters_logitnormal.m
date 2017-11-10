function [parameter_out, changes] = update_parameters_logitnormal(...
    parameter,loglklh, logprior, logvariational,...
    dqdp_logit, dqdalpha, dqdbeta, dqdalpha_gain, dqdbeta_gain,...
    iteration,eta,eta_max,varargin)


if ~isempty(varargin) && ~isempty(varargin{1})
   spike_indicator= varargin{1};
else
    spike_indicator= false;
end


n_cell=size(dqdalpha,1);
S=size(dqdalpha,2);

f_p_logit=zeros(n_cell,S);
f_alpha=zeros(n_cell,S);f_beta=zeros(n_cell,S);
f_alpha_gain=zeros(n_cell,S);f_beta_gain=zeros(n_cell,S);
h_p_logit=zeros(n_cell,S);
h_alpha=zeros(n_cell,S);h_beta=zeros(n_cell,S);
h_alpha_gain=zeros(n_cell,S);h_beta_gain=zeros(n_cell,S);


parameter_out=parameter;
    changes=0;
for i_cell = 1:n_cell
    sum_of_logs = loglklh(i_cell,:)+logprior(i_cell,:)-logvariational(i_cell,:);
    
    % need to wrap up, either vectorize or as one function:
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
    
    
    a_constant = (spike_indicator*quick_cov(f_p_logit(i_cell,:),h_p_logit(i_cell,:))+quick_cov(f_alpha(i_cell,:),h_alpha(i_cell,:))+...
        quick_cov(f_beta(i_cell,:),h_beta(i_cell,:))+...
        quick_cov(f_alpha_gain(i_cell,:),h_alpha_gain(i_cell,:))+...
        quick_cov(f_beta_gain(i_cell,:),h_beta_gain(i_cell,:))...
        )/( spike_indicator*quick_cov(h_p_logit(i_cell,:),h_p_logit(i_cell,:))+ ...
        quick_cov(h_alpha(i_cell,:),h_alpha(i_cell,:))+...
        quick_cov(h_beta(i_cell,:),h_beta(i_cell,:))+...
        quick_cov(h_alpha_gain(i_cell,:),h_alpha_gain(i_cell,:))+...
        quick_cov(h_beta_gain(i_cell,:),h_beta_gain(i_cell,:))...
        +0.01);
    
    step_size = (eta/iteration);
    grad_logit=spike_indicator*step_size*mean(f_p_logit(i_cell,:)-a_constant*h_p_logit(i_cell,:));
    grad_alpha = step_size*mean(f_alpha(i_cell,:)-a_constant*h_alpha(i_cell,:));
    grad_beta = step_size*mean(f_beta(i_cell,:)-a_constant*h_beta(i_cell,:));
    grad_alpha_gain=step_size*mean(f_alpha_gain(i_cell,:)-a_constant*h_alpha_gain(i_cell,:));
    grad_beta_gain=step_size*mean(f_beta_gain(i_cell,:)-a_constant*h_beta_gain(i_cell,:));
    grad_max = max(abs([grad_logit grad_alpha grad_beta grad_alpha_gain grad_beta_gain]));
    
    if grad_max > eta_max
        grad_scale= grad_max/eta_max;
        grad_logit= grad_logit/grad_scale;
        grad_alpha = grad_alpha/grad_scale;
        grad_beta = grad_beta/grad_scale;
        grad_alpha_gain=grad_alpha_gain/grad_scale;
        grad_beta_gain=grad_beta_gain/grad_scale;
    end
    
    
    parameter_out(i_cell).p_logit=parameter(i_cell).p_logit+spike_indicator*grad_logit;
    parameter_out(i_cell).alpha=parameter(i_cell).alpha+grad_alpha;
    parameter_out(i_cell).beta=parameter(i_cell).beta+grad_beta;
    parameter_out(i_cell).alpha_gain=parameter(i_cell).alpha_gain+grad_alpha_gain;
    parameter_out(i_cell).beta_gain=parameter(i_cell).beta_gain+grad_beta_gain;
    changes= changes+  abs(grad_logit)+abs(grad_alpha)+abs(grad_beta)+abs(grad_alpha_gain)+...
        abs(grad_beta_gain);
    
    
end

    changes= changes/(sum(abs([parameter(:).p_logit]))+...
        sum(abs([parameter(:).alpha]))+sum(abs([parameter(:).beta]))+...
        sum(abs([parameter(:).alpha_gain]))+sum(abs([parameter(:).beta_gain])));
    
end
