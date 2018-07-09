function [parameter_history, elbo_rec,loglklh_rec] = fit_VI_1D(...
    trials,neurons, background_rate, ...
    variational_params,prior_params,...
    inference_params,prior_info,GP_samples)
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
boundary=80;
%%
spike_curves=prior_info.induced_intensity;
parameter_current=variational_params;

iteration = 1;loglklh_rec=[];
%%
tic;
% time_rec=zeros(10,1);
while (change_history(iteration) > epsilon && iteration<maxit)
    %%
    tstart=toc;
    parameter_history(iteration,:)=parameter_current;
    iteration=iteration+1;
    adjusted_location=zeros(n_cell,1);
    loglklh=zeros(S,1);logprior=zeros(S,1);logvariational=zeros(S,1);
    %% Calculate the stim size using samples from GP and the shifts
    vsam=cell([S 1]);rsam=cell([S 1]);
    for s=1:S
        [variational_samples,raw_samples] = draw_samples_from_var_dist(parameter_current);
        vsam{s}=variational_samples;rsam{s}=raw_samples;
        logprior(s)=get_logdistribution(variational_samples,raw_samples,prior_params);
        logvariational(s)=get_logdistribution(variational_samples,raw_samples,parameter_current);
        [loglklh(s)] = update_likelihood(trials, variational_samples,...
            GP_samples, background_rate,lklh_func,spike_curves,neurons,boundary);
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
    new_gradient=sum_gradient(gradients,eta,eta_max,iteration);
    [parameter_current, change_history(iteration)]=incorporate_gradient(parameter_current, new_gradient);
    loglklh_rec(iteration)=mean(mean(loglklh));
    elbo_rec(iteration)=mean(logprior+loglklh-logvariational);
    tdiff=toc-tstart;
    fprintf('Iteration %d; change %d; ELBO %d; time %d; \n',iteration,change_history(iteration),elbo_rec(iteration),tdiff)
end
fprintf('VI fitted after %d iteration;\n',iteration)

