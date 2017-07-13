function out = derinvlink_sig(resp)

exp_part = exp(-3*(resp+1));
out = (1-0.01)*3*exp_part/((1+exp_part)^2);