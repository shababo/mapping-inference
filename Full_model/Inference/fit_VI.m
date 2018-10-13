function [parameter_history, elbo_rec] = fit_VI(...
    trials,neurons, background_rate, ...
    variational_params,prior_params,...
    inference_params,prior_info)
%%
if isfield(variational_params(1),'PR')
par_gradients=true;
else 
par_gradients=false;
    
end
    
    
n_cell=length(neurons);
n_trial = length(trials);

epsilon=inference_params.convergence_threshold;
S=inference_params.MCsamples_for_gradient;
eta=inference_params.step_size;
eta_max=inference_params.step_size_max;
eta_threshold=inference_params.eta_threshold;
maxit=inference_params.maxit;
lklh_func=inference_params.likelihood;
change_history=epsilon+1;
boundary_params = prior_info.prior_parameters.boundary_params;
GP_params=prior_info.prior_parameters.GP_params;
spike_curves=prior_info.induced_intensity;
parameter_current=variational_params;

if isfield(inference_params, 'event_range')
   for i_trial = 1:length(trials)
      trials(i_trial).event_times=trials(i_trial).event_times( trials(i_trial).event_times < inference_params.event_range(2) &...
          trials(i_trial).event_times> inference_params.event_range(1));
   end
end
iteration = 1;loglklh_rec=[];

% Calculate a grid for standard normal r.v. for quick CDF calculation:
pre_density=struct;
grid.bound=4;grid.gap=0.1;
pre_density.grid=grid;
pre_density.normal_grid = -grid.bound:grid.gap:grid.bound;
for i = 1:length(pre_density.normal_grid)
    pre_density.cdf_grid(i) = normcdf(pre_density.normal_grid(i),0,1);
    pre_density.pdf_grid(i)=normpdf(pre_density.normal_grid(i),0,1);
end
% S=100;
% clear('parameter_history')
% clear('gradients')
%%
if par_gradients
loglklh_PR=zeros(S,1);logprior_PR=zeros(S,1);logvariational_PR=zeros(S,1);

    [variational_samples,raw_samples] = draw_samples_from_var_dist(parameter_current);
       
variational_mean =variational_samples;
fldnames= fieldnames(parameter_current);
raw_mean = variational_samples;
for i_neuron =1:n_cell
    for j= 1:length(fldnames)
        raw_mean(i_neuron).(fldnames{j})=variational_params(i_neuron).(fldnames{j}).mean;
%         if ~strcmp(fldnames{j},'shapes')
            this_raw = raw_mean(i_neuron).(fldnames{j});
            lb= variational_params(i_neuron).(fldnames{j}).bounds.low;
            ub= variational_params(i_neuron).(fldnames{j}).bounds.up;
            variational_mean(i_neuron).(fldnames{j})= exp(this_raw)./(1+exp(this_raw)).*(ub-lb)+lb;
            
%         end
    end
end
end 
%%
tic;
% ts=zeros(10,1);
while (change_history(iteration) > epsilon && iteration<maxit)
    %%
    tstart=toc;
    parameter_history(iteration,:)=parameter_current;
    iteration=iteration+1;
    adjusted_location=zeros(n_cell,3);
    loglklh=zeros(S,1);logprior=zeros(S,1);logvariational=zeros(S,1);
    %% Calculate the stim size using samples from GP and the shifts
    if ~par_gradients
        vsam=cell([S 1]);rsam=cell([S 1]);
        for s=1:S
            %         t1=toc;
            [variational_samples,raw_samples] = draw_samples_from_var_dist(parameter_current);
            
            vsam{s}=variational_samples;rsam{s}=raw_samples;
            logprior(s)=get_logdistribution(variational_samples,raw_samples,prior_params);
            logvariational(s)=get_logdistribution(variational_samples,raw_samples,parameter_current);
            [loglklh(s)] = update_likelihood(trials, variational_samples,parameter_current,...
                background_rate,lklh_func,spike_curves,neurons,prior_info,inference_params,pre_density);
            
            lklhweight=logprior(s)+loglklh(s)-logvariational(s);
            this_gradient=get_variational_gradient(variational_samples,raw_samples, parameter_current);
            this_gradient=get_variate_control(lklhweight,this_gradient);
            if exist('gradients')
                gradients(s,:) = this_gradient;
            else
                gradients(s,:)=this_gradient;
            end
        end
    else
        vsam=cell([S 1]);rsam=cell([S 1]);
        for s=1:S
            %         t1=toc;
            [variational_samples,raw_samples] = draw_samples_from_var_dist(parameter_current);
            %% Calculate the gradients for PR and other parameters separately
            vsam{s}=variational_samples;rsam{s}=raw_samples;
            % for PR:
            var_sample=  variational_mean;
            r_sample = raw_mean;
            for i_neuron =1:n_cell
                var_sample(i_neuron).PR=variational_samples(i_neuron).PR;
                r_sample(i_neuron).PR=raw_samples(i_neuron).PR;
            end
            logprior_PR(s)=get_logdistribution(var_sample,r_sample,prior_params);
            logvariational_PR(s)=get_logdistribution(var_sample,r_sample,parameter_current);
            [loglklh_PR(s)] = update_likelihood(trials, var_sample,parameter_current,...
                background_rate,lklh_func,spike_curves,neurons,prior_info,inference_params,pre_density);
            lklhweight=logprior_PR(s)+loglklh_PR(s)-logvariational_PR(s);
            this_gradient_PR=get_variational_gradient(var_sample,r_sample, parameter_current);
            this_gradient_PR=get_variate_control(lklhweight,this_gradient_PR);
            % For the rest of parameters 
            for i_neuron =1:n_cell
                variational_samples(i_neuron).PR=variational_mean(i_neuron).PR;
                raw_samples(i_neuron).PR=raw_mean(i_neuron).PR;
            end
          
            logprior(s)=get_logdistribution(variational_samples,raw_samples,prior_params);
            logvariational(s)=get_logdistribution(variational_samples,raw_samples,parameter_current);
            [loglklh(s)] = update_likelihood(trials, variational_samples,parameter_current,...
                background_rate,lklh_func,spike_curves,neurons,prior_info,inference_params,pre_density);
            lklhweight=logprior(s)+loglklh(s)-logvariational(s);
            this_gradient=get_variational_gradient(variational_samples,raw_samples, parameter_current);
            this_gradient=get_variate_control(lklhweight,this_gradient);
            
            % merge gradients: 
            for i_neuron = 1:n_cell
            this_gradient(i_neuron).PR=this_gradient_PR(i_neuron).PR;
            end
            
            if exist('gradients')
                gradients(s,:) = this_gradient;
            else
                gradients(s,:)=this_gradient;
            end
        end
    end
    %% 
    
%     gains=[];PRs=[];sigmafs=[];meanfs=[];
%     for i = 1:S
%     gains(i) = vsam{i}(1).gain; 
% %     PRs(i) = vsam{i}.PR; 
% %     sigmafs(i)=gradients(i).PR.sigma_f;
% %     meanfs(i)=gradients(i).PR.mean_f;
% %     lklhweight(i)=logprior(i)+loglklh(i)-logvariational(i);
%     end
%     figure(1)
%     scatter(gains,loglklh)
% %     figure(2)
% %     scatter(meanfs,loglklh) 
%     figure(2)
%     scatter(PRs,loglklh)
%     figure(3)
    
%     scatter(PRs,lklhweight)
%      figure(3)
%     scatter(PRs,meanfs)
%     figure(4)
%     scatter(PRs,sigmafs)
%     mean(sigmafs)
% mean(meanfs)

% traces=[];traces2=[];
% for i = 1:(iteration-1)
% traces(i)=parameter_history(i).PR.mean;
% traces2(i)=parameter_history(i).PR.log_sigma;
% end
% figure(1)
% plot(traces)
% 
% figure(2)
% plot(traces2)
% figure(3)
% plot(elbo_rec(3:iteration))
%%
    new_gradient=sum_gradient(gradients,eta,eta_max,iteration,eta_threshold);
    [parameter_current, change_history(iteration)]=incorporate_gradient(parameter_current, new_gradient);
    elbo_rec(iteration)=mean(logprior+loglklh-logvariational);
    
    tend=toc;
    tdiff=tend-tstart;
    fprintf('Iteration %d; change %d; time %d; \n',iteration,change_history(iteration),tdiff)
    % ELBO %d;PR %d; ,elbo_rec(iteration),parameter_current(1,1).PR.mean)
%     fprintf('PR %d; \n',parameter_current.PR.mean)
end
fprintf('VI fitted after %d iteration;\n',iteration)
