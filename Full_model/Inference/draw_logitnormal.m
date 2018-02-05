function [this_param_sample,this_raw_sample]= draw_logitnormal(this_param, S)
    alpha=this_param.params(1);
    beta=this_param.params(2);
    this_raw_sample=normrnd(alpha,exp(beta),[1 S]);
    this_param_sample=exp(this_raw_sample);
    this_param_sample=this_param_sample./(1+this_param_sample);

    
