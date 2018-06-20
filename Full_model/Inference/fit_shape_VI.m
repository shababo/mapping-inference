function [parameter_history, loglklh_rec,elbo_rec] = fit_shape_VI(...
    neurons,variational_params,prior_params,inference_params,prior_opt)
%%
epsilon=inference_params.convergence_threshold;
S=inference_params.MCsamples_for_gradient;
eta=inference_params.step_size;
eta_max=inference_params.step_size_max;
maxit=inference_params.maxit;
mean_func=inference_params.mean_func;
n_cell = length(neurons);
% initialize storages
logvariational=zeros(S,1);
logprior=zeros(S,1);loglklh=zeros(S,1);

corrected_grid=cell([S n_cell]);
all_y=[neurons(:).scaled_current]';

parameter_history=variational_params;

parameter_current=variational_params;

iteration = 1;loglklh_rec=[];
elbo_rec=0;
change_history(iteration) = epsilon+1;
% lklh_func=inference_params.likelihood;

%%
% tic;
% time_rec=zeros(10,1);
while (change_history(iteration) > epsilon && iteration<maxit)
    %%
    parameter_history(iteration,:)=parameter_current;
    iteration=iteration+1;
    clear('gradients')
    for s= 1:S
        variational_samples = draw_samples_from_var_dist(parameter_current);
        % Calculate the new grid
        
        for i_cell = 1:n_cell
            corrected_grid{s,i_cell}=[neurons(i_cell).stim_grid - variational_samples(i_cell).shift];
        end
       
        % Calculate the log pdf of prior  dist.
        if prior_opt
        logprior(s)=get_logdistribution(variational_samples,prior_params);
        else
        logprior(s)=0;
        end
        % Calculate the log pdf of variational dist.
        logvariational(s)=get_logdistribution(variational_samples,parameter_current);
        
        loglklh(s)=0;
        for i_cell = 1:n_cell
            X=corrected_grid{s,i_cell}';
            Y=neurons(i_cell).scaled_current';
            nsq=sum(X.^2,2);
            K=bsxfun(@plus,nsq,nsq');
            K=bsxfun(@minus,K,(2*X)*X.');
            K=variational_samples(1).GP_sigma*exp(-K/variational_samples(1).GP_tau);
            K=K+ diag(ones(length(corrected_grid{s,i_cell}),1))*variational_samples(i_cell).current_sigma;
            loglklh(s)=loglklh(s)+log(mvnpdf(Y,mean_func(X),K));
        end
        
        
        if isinf(loglklh(s))
           loglklh(s)= -1e5; 
        end
        
        lklhweight = logprior(s)+loglklh(s)-logvariational(s);
        % 
        this_gradient=get_variational_gradient(variational_samples,parameter_current);
        
        %obtain the f and h functions: 
        this_gradient=get_variate_control(lklhweight,this_gradient);
        
        if exist('gradients')
           gradients(s,:) = this_gradient;
        else
           gradients(s,:)=this_gradient;
        end
    end
    
    new_gradient=sum_gradient(gradients,eta,eta_max,iteration);
    
    [parameter_current, change_history(iteration)]=incorporate_gradient(parameter_current, new_gradient);
    
    
    loglklh_rec(iteration)=mean(mean(loglklh));
    elbo_rec(iteration)=mean(logprior+loglklh-logvariational);
   
    fprintf('Iteration %d; change %d; ELBO %d \n',iteration,change_history(iteration),elbo_rec(iteration))
end
fprintf('VI fitted after %d iteration;\n',iteration)
%%
% exp(parameter_current(1).GP_tau.mean)

% exp(parameter_current(1).GP_sigma.mean)

% exp(parameter_current(1).current_sigma.mean)

% plot(loglklh_rec)
