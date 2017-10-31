function [logit_sample, sample] = draw_spiked_logitnormal(alpha,beta,bounds,varargin)
n_cell=length(alpha);
if ~isempty(varargin) && ~isempty(varargin{1})
   zero_prob=exp(varargin{1})./(exp(varargin{1})+1); 
   zero_prob=zero_prob';
else
     zero_prob=zeros(n_cell,1);
end


if length(varargin)>1 && ~isempty(varargin{2})
  spike_indicator=varargin{2};
else
  spike_indicator=false;
end

if ~spike_indicator
    zero_prob=zeros(n_cell,1); 
end

beta=exp(beta);
spike_sample = rand(n_cell,1) > (zero_prob);
logit_sample=normrnd(alpha',beta',[n_cell 1]);
slab_sample = exp(logit_sample)./(1+exp(logit_sample))*(bounds.up-bounds.low)+bounds.low;
sample = slab_sample.*spike_sample;
end




