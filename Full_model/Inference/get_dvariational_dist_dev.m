function [dlogit, dalpha,dbeta] = get_dvariational_dist(sample,logit_sample,alpha,beta,constants,varargin)
n_cell=length(alpha);
if ~isempty(varargin) &&  ~isempty(varargin{1})
    p_logit=varargin{1}';
else
    p_logit=zeros(n_cell,1);
end

if length(varargin)>1 &&  ~isempty(varargin{2})
     spike_indicator=varargin{2};
else
  spike_indicator=false;
end

if spike_indicator
dlogit= (sample==0)./(1+exp(p_logit))-...
    (sample>0).*exp(p_logit)./(1+exp(p_logit));
else
dlogit= 0*beta;   
end


beta=exp(beta);
dalpha= (-constants.mu_by_sigma2+logit_sample.*constants.sigma_inv2);
% dalpha(sample==0)=0;
dbeta= (-1+constants.sigma_inv2.*(logit_sample-alpha').^2);
% dbeta(sample==0)=0;

end