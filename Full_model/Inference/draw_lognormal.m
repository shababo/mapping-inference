function [this_param_sample,this_raw_sample]= draw_lognormal(this_param, S)
    alpha=this_param.params(1);
    beta=this_param.params(2);
    this_raw_sample=normrnd(alpha,exp(beta),[1 S]);
    this_param_sample=exp(this_raw_sample);

    
