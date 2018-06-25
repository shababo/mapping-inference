function [parameter_history, loglklh_rec,elbo_rec] = fit_shape_VI(...
    neurons,variational_params,prior_params,inference_params)
%%
epsilon=inference_params.convergence_threshold;
S=inference_params.MCsamples_for_gradient;
eta=inference_params.step_size;
eta_max=inference_params.step_size_max;
maxit=inference_params.maxit;
mean_func=inference_params.mean_func;
var_func=inference_params.var_func;
n_cell = length(neurons);

prior_opt =inference_params.prior_opt;
%     noise_sigma =[neurons(:).noise_sigma];

% initialize storages
eig_epsilon=inference_params.eig_epsilon;
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
            corrected_grid{s,i_cell}=[neurons(i_cell).adjusted_grid - variational_samples(i_cell).shift];
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
            mean_func=inference_params.mean_func;
            
            % What is our variance estimator if we were to use the unscaled
            % version?
            sigma_grid = sqrt(max(eig_epsilon,var_func.func(X,var_func.params)));
            
            Y=neurons(i_cell).scaled_current';
            white_sigma = neurons(i_cell).noise_sigma;
            Ymean = mean_func.func(X,mean_func.params);
            if inference_params.fit_gain
                
                Y=neurons(i_cell).raw_current';
                sigma_grid = variational_samples(i_cell).gain*sigma_grid;
                white_sigma =  variational_samples(i_cell).gain*white_sigma;
                Ymean=Ymean*variational_samples(i_cell).gain;
            end
            nsq=sum(X.^2,2);
            K=bsxfun(@plus,nsq,nsq');
            K=bsxfun(@minus,K,(2*X)*X.');
            K=exp(-K/variational_samples(1).GP_tau);
            sigma_mat = sigma_grid*ones(1,length(X));
            Kcov=sigma_mat.*K.*sigma_mat';
            Kcov=(Kcov+Kcov')/2;
            
            Kmarcov=Kcov+ diag(white_sigma.^2);
            
            loglklh(s)=loglklh(s)+log(mvnpdf(Y,Ymean,Kmarcov));
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
