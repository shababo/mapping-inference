function [out]=expit(x,lwb,upb)
out=[exp(x)./(1+exp(x))].*(upb-lwb)+lwb;
