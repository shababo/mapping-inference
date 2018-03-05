function [delay_params] = ...
    get_Delay_pilot(mpp,stim_size,spike_curves,varargin)
%%
% FOr eahc cell:


if ~isempty(varargin) && ~isempty(varargin{1})
    
else
   
end
%%
[power_levels,ia,ic] = unique(stim_size);
mean_vec=zeros(length(power_levels),1);
var_vec=zeros(length(power_levels),1);
for i=1:length(power_levels)
    trial_index=find(ic==i);
    mean_vec(i)=mean([mpp(trial_index).event_times]);
    var_vec(i)=var([mpp(trial_index).event_times]);
end


ni_mean=isnan(mean_vec);
ydata=mean_vec(~ni_mean);
ydata2=var_vec(~ni_mean);
xdata=power_levels(~ni_mean)'*1e3;
convo_mean = @(x,xdata) x(1) + spike_curves.specs.F_mean(spike_curves.mean_param,x(2).*xdata);
convo_var=@(x,xdata) x(3) +  ...
    spike_curves.specs.F_sd(spike_curves.sd_param,x(2).*xdata).^2;
x0 = [20 0.05 10];
Fsumsquares = @(x)sum((convo_mean(x,xdata) - ydata).^2 + (convo_var(x,xdata) - ydata2).^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton');
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0 0 0];
ub = [];
% fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon) ,opts
[convo_param,ressquared,eflag,outputu] = fmincon(Fsumsquares,x0,A,b,Aeq,beq,lb,ub);
%convo_param(1) delay mean
% convo_param(2) gain
    convo_param
    %% Obtain a smooth curve the for errors in fitting a mean curve
    % and Obtain a crude curve for SD (within condition variability)
    
     %[~, stim_index]=min(abs(effective_stim - spike_curves.current));
       
     avai_index=find(var_vec>0);
     xdata=convo_param(2)*power_levels(avai_index);
     
     act_curr_index = zeros(length(xdata),1);
     var_explained = zeros(length(xdata),1);
     
     for i =1:length(xdata)
        [~, act_curr_index(i)]=min(abs(xdata(i) - spike_curves.current));
         var_explained(i)=spike_curves.sd(act_curr_index(i))^2+spike_curves.dev(act_curr_index(i))^2;
     end
     
     
    %% Estimate the variance of convo_param based on the fitted means:
    delay_params=struct;
    
    delay_params.var=convo_param(3);
    delay_params.mean=convo_param(1);
    delay_params.gain=convo_param(2);
    
    delay_params.mpp=mpp;
    delay_params.stim=stim_size;
    