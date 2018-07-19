function [parameter_history, elbo_rec] = fit_VI(...
    trials,neurons, background_rate, ...
    variational_params,prior_params,...
    inference_params,prior_info)
%%
n_cell=length(neurons);
n_trial = length(trials);

epsilon=inference_params.convergence_threshold;
S=inference_params.MCsamples_for_gradient;
eta=inference_params.step_size;
eta_max=inference_params.step_size_max;
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
% S=50;
% clear('parameter_history')
% clear('gradients')
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
    vsam=cell([S 1]);rsam=cell([S 1]);
    for s=1:S
%         t1=toc;
        [variational_samples,raw_samples] = draw_samples_from_var_dist(parameter_current);
%         if only_PR
%             variational_samples.gain=neurons.truth.optical_gain;
%             variational_samples.delay_mu=neurons.truth.delay_mean;
%             variational_samples.delay_sigma=sqrt(neurons.truth.delay_var);
%             variational_samples.shapes=0.747;
%         else
% %             variational_samples.PR=0.5; % fixed PR
%         end
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
    %%
%     gains=[];PRs=[];sigmafs=[];meanfs=[];
%     for i = 1:S
%     gains(i) = vsam{i}.gain; 
%     PRs(i) = vsam{i}.PR; 
%     sigmafs(i)=gradients(i).PR.sigma_f;
%     meanfs(i)=gradients(i).PR.mean_f;
%     lklhweight(i)=logprior(i)+loglklh(i)-logvariational(i);
%     end
%     figure(1)
%     scatter(gains,loglklh)
%     figure(2)
%     scatter(meanfs,loglklh) 
%     figure(2)
%     scatter(PRs,loglklh)
%     figure(3)
%     
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
    new_gradient=sum_gradient(gradients,eta,eta_max,iteration);
    [parameter_current, change_history(iteration)]=incorporate_gradient(parameter_current, new_gradient);
    elbo_rec(iteration)=mean(logprior+loglklh-logvariational);
    tend=toc;
    tdiff=tend-tstart;
    fprintf('Iteration %d; change %d; ELBO %d; time %d; \n',iteration,change_history(iteration),elbo_rec(iteration),tdiff)
%     fprintf('PR %d; \n',parameter_current.PR.mean)
end
fprintf('VI fitted after %d iteration;\n',iteration)
