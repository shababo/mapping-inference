function [output] = draw_likelihood_samples(...
    params,trials,neurons, background_rate, ...
    variational_params,...
    inference_params,prior_info)
%%
draw_truth=false;
if isfield(params,'truth')
    draw_truth=true;
end

%%
S=500;
eta=1;
eta_max=1;
eta_threshold=100;
iteration=1;
lklh_func=inference_params.likelihood;
boundary_params = prior_info.prior_parameters.boundary_params;
GP_params=prior_info.prior_parameters.GP_params;
spike_curves=prior_info.induced_intensity;

if isfield(inference_params, 'event_range')
   for i_trial = 1:length(trials)
      trials(i_trial).event_times=trials(i_trial).event_times( trials(i_trial).event_times < inference_params.event_range(2) &...
          trials(i_trial).event_times> inference_params.event_range(1));
   end
end

% Calculate a grid for standard normal r.v. for quick CDF calculation:
pre_density=struct;
grid.bound=4;grid.gap=0.1;
pre_density.grid=grid;
pre_density.normal_grid = -grid.bound:grid.gap:grid.bound;
for i = 1:length(pre_density.normal_grid)
    pre_density.cdf_grid(i) = normcdf(pre_density.normal_grid(i),0,1);
    pre_density.pdf_grid(i)=normpdf(pre_density.normal_grid(i),0,1);
end
%% Generate samples based on given information 
variational_samples=cell([S 1]);raw_samples=cell([S 1]);
for s=1:S
    [variational_samples{s},raw_samples{s}] = draw_samples_from_var_dist(variational_params);
end
%% Calculate loglikelihoods using these samples 
loglklh=zeros(S,1);logprior=zeros(S,1);logvariational=zeros(S,1);
for s=1:S
    logprior(s)=0;
    logvariational(s)=get_logdistribution(variational_samples{s},raw_samples{s},variational_params);
    [loglklh(s)] = update_likelihood(trials, variational_samples{s},variational_params,...
        background_rate,lklh_func,spike_curves,neurons,prior_info,inference_params,pre_density);
    lklhweight=logprior(s)+loglklh(s)-logvariational(s);
    this_gradient=get_variational_gradient(variational_samples{s},raw_samples{s}, variational_params);
    this_gradient=get_variate_control(lklhweight,this_gradient);
    if exist('gradients')
        gradients(s,:) = this_gradient;
    else
        gradients(s,:)=this_gradient;
    end
end
  new_gradient=sum_gradient(gradients,eta,eta_max,iteration,eta_threshold);
  
%% Visualize the log-likelihood as functions of samples 

% first, turn the cell array into vectors for each field and each neuron
n_neurons = length(variational_samples{1});
clear('neuron_samples');
neuron_samples(n_neurons)=struct;

fldnames=fieldnames(variational_params);
for i_neuron= 1:n_neurons
    for j= 1:length(fldnames)
        if ~strcmp(fldnames{j},'shapes')
            this_sample=zeros(S,1);
            for s=1:S
                this_sample(s)=variational_samples{s}(i_neuron).(fldnames{j});
                
            end
            neuron_samples(i_neuron).(fldnames{j})=this_sample;
        end
    end
end
%% Now draw scatter plots between loglklh and parameters 
if params.plot.do
 figure(1)            
for i_neuron= 1:n_neurons
    for j= 1:length(fldnames)
        if ~strcmp(fldnames{j},'shapes')
            this_sample=neuron_samples(i_neuron).(fldnames{j});
            i_plot = (i_neuron-1)*(length(fldnames)-1)+j;
            subplot(n_neurons,length(fldnames)-1,i_plot);
            scatter(this_sample,loglklh,'MarkerFaceColor','b')
            hold on;
            if draw_truth
                true_value = params.truth.neurons(i_neuron).truth.(fldnames{j})*ones(2,1);
                line(true_value, [min(loglklh) max(loglklh)],'color','r','LineStyle','--')
                hold on;
                title_string=['Neuron ' num2str(i_neuron) '; PR = ' num2str(params.truth.neurons(i_neuron).truth.PR ) ];
            end
            xlabel(fldnames{j})
            ylabel('Log-Likelihood')
            title(title_string)
            % Draw the current mean and the direction of gradients:
            this_param_raw=variational_params(i_neuron).(fldnames{j}).mean;
            lb= variational_params(i_neuron).(fldnames{j}).bounds.low;
            ub= variational_params(i_neuron).(fldnames{j}).bounds.up;
            
            this_param= exp(this_param_raw)/(1+exp(this_param_raw))*(ub-lb)+lb;
            
            this_grad = new_gradient(i_neuron).(fldnames{j}).mean;
            new_param_raw=this_param_raw+ this_grad;
            new_param =  exp(new_param_raw)/(1+exp(new_param_raw))*(ub-lb)+lb;
               
            norm_diff = 0.15*range(this_sample)*(new_param-this_param)/abs((new_param-this_param));
                line(this_param*ones(2,1), [min(loglklh) max(loglklh)],'color','g','LineStyle','-')
                hold on;
                % show the direction of gradients:
                quiver(this_param,mean(loglklh),norm_diff,0,'AutoScale','off', 'Color','g','LineWidth',2,'Marker','o','MarkerFaceColor','g'); 
%                 annotation('textarrow',[this_param  mean(loglklh)],[this_param+norm_diff  mean(loglklh)],...
%                     'Color','b','LineWidth',2);
                
        end
    end
end
end 
%%
output=struct;
output.gradients=gradients;
output.new_gradient=new_gradient;
output.variational_samples=variational_samples;
output.raw_samples=raw_samples;
output.loglklh=loglklh;

%% New strategies for estimationg the loglikelihood & gradients:
loglklh_PR=zeros(S,1);logprior_PR=zeros(S,1);logvariational_PR=zeros(S,1);
loglklh_others=zeros(S,1);logprior_others=zeros(S,1);logvariational_others=zeros(S,1);
variational_mean =variational_samples{1};

raw_mean = variational_samples{1};

for i_neuron =1:n_neurons
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
% 
%%
 for j= 1:length(fldnames)
        if strcmp(fldnames{j},'PR')
            for s=1:S
            var_sample=  variational_mean;
            r_sample = raw_mean;
            for i_neuron =1:n_neurons
                var_sample(i_neuron).PR=variational_samples{s}(i_neuron).PR;
                r_sample(i_neuron).PR=raw_samples{s}(i_neuron).PR;
            end
            logprior_PR(s)=0;
            logvariational_PR(s)=get_logdistribution(var_sample,r_sample,variational_params);
            [loglklh_PR(s)] = update_likelihood(trials, var_sample,variational_params,...
                background_rate,lklh_func,spike_curves,neurons,prior_info,inference_params,pre_density);
            lklhweight=logprior_PR(s)+loglklh_PR(s)-logvariational_PR(s);
            this_gradient_PR=get_variational_gradient(var_sample,r_sample, variational_params);
            this_gradient_PR=get_variate_control(lklhweight,this_gradient_PR);
            if exist('gradients_PR')
                gradients_PR(s,:) = this_gradient_PR;
            else
                gradients_PR(s,:)=this_gradient_PR;
            end
            end
        else
             for s=1:S
            var_sample=  variational_samples{s};
            r_sample = raw_samples{s};
            for i_neuron =1:n_neurons
                var_sample(i_neuron).PR=variational_mean(i_neuron).PR;
                r_sample(i_neuron).PR=raw_mean(i_neuron).PR;
            end
            logprior_others(s)=0;
            logvariational_others(s)=get_logdistribution(var_sample,r_sample,variational_params);
            [loglklh_others(s)] = update_likelihood(trials, var_sample,variational_params,...
                background_rate,lklh_func,spike_curves,neurons,prior_info,inference_params,pre_density);
            lklhweight=logprior_others(s)+loglklh_others(s)-logvariational_others(s);
            this_gradient_others=get_variational_gradient(var_sample,r_sample, variational_params);
            this_gradient_others=get_variate_control(lklhweight,this_gradient_others);
            if exist('gradients_others')
                gradients_others(s,:) = this_gradient_others;
            else
                gradients_others(s,:)=this_gradient_others;
            end
            end
            
        end
 end
 
  new_gradient_PR=sum_gradient(gradients_PR,eta,eta_max,iteration,eta_threshold);
  new_gradient_others=sum_gradient(gradients_others,eta,eta_max,iteration,eta_threshold);
  

%% Visualize the new gradients
if params.plot.do
figure(2)         
for i_neuron= 1:n_neurons
    for j= 1:length(fldnames)
        if strcmp(fldnames{j},'PR')
            this_sample=neuron_samples(i_neuron).(fldnames{j});
            i_plot = (i_neuron-1)*(length(fldnames)-1)+j;
            subplot(n_neurons,length(fldnames)-1,i_plot);
            scatter(this_sample,loglklh_PR,'MarkerFaceColor','b')
            hold on;
            if draw_truth
                true_value = params.truth.neurons(i_neuron).truth.(fldnames{j})*ones(2,1);
                line(true_value, [min(loglklh_PR) max(loglklh_PR)],'color','r','LineStyle','--')
                hold on;
                title_string=['Neuron ' num2str(i_neuron) '; PR = ' num2str(params.truth.neurons(i_neuron).truth.PR ) ];
            end
            xlabel(fldnames{j})
            ylabel('Log-Likelihood')
            title(title_string)
            % Draw the current mean and the direction of gradients:
            this_param_raw=variational_params(i_neuron).(fldnames{j}).mean;
            lb= variational_params(i_neuron).(fldnames{j}).bounds.low;
            ub= variational_params(i_neuron).(fldnames{j}).bounds.up;
            
            this_param= exp(this_param_raw)/(1+exp(this_param_raw))*(ub-lb)+lb;
            
            this_grad = new_gradient_PR(i_neuron).(fldnames{j}).mean;
            new_param_raw=this_param_raw+ this_grad;
            new_param =  exp(new_param_raw)/(1+exp(new_param_raw))*(ub-lb)+lb;
               
            norm_diff = 0.15*range(this_sample)*(new_param-this_param)/abs((new_param-this_param));
                line(this_param*ones(2,1), [min(loglklh_PR) max(loglklh_PR)],'color','g','LineStyle','-')
                hold on;
                % show the direction of gradients:
                quiver(this_param,mean(loglklh_PR),norm_diff,0,'AutoScale','off', 'Color','g','LineWidth',2,'Marker','o','MarkerFaceColor','g'); 
               
        elseif ~strcmp(fldnames{j},'shapes')
            
            this_sample=neuron_samples(i_neuron).(fldnames{j});
            i_plot = (i_neuron-1)*(length(fldnames)-1)+j;
            subplot(n_neurons,length(fldnames)-1,i_plot);
            scatter(this_sample,loglklh_others,'MarkerFaceColor','b')
            hold on;
            if draw_truth
                true_value = params.truth.neurons(i_neuron).truth.(fldnames{j})*ones(2,1);
                line(true_value, [min(loglklh_others) max(loglklh_others)],'color','r','LineStyle','--')
                hold on;
                title_string=['Neuron ' num2str(i_neuron) '; PR = ' num2str(params.truth.neurons(i_neuron).truth.PR ) ];
            end
            xlabel(fldnames{j})
            ylabel('Log-Likelihood')
            title(title_string)
            % Draw the current mean and the direction of gradients:
            this_param_raw=variational_params(i_neuron).(fldnames{j}).mean;
            lb= variational_params(i_neuron).(fldnames{j}).bounds.low;
            ub= variational_params(i_neuron).(fldnames{j}).bounds.up;
            
            this_param= exp(this_param_raw)/(1+exp(this_param_raw))*(ub-lb)+lb;
            
            this_grad = new_gradient_others(i_neuron).(fldnames{j}).mean;
            new_param_raw=this_param_raw+ this_grad;
            new_param =  exp(new_param_raw)/(1+exp(new_param_raw))*(ub-lb)+lb;
               
            norm_diff = 0.15*range(this_sample)*(new_param-this_param)/abs((new_param-this_param));
                line(this_param*ones(2,1), [min(loglklh_others) max(loglklh_others)],'color','g','LineStyle','-')
                hold on;
                % show the direction of gradients:
                quiver(this_param,mean(loglklh_others),norm_diff,0,'AutoScale','off', 'Color','g','LineWidth',2,'Marker','o','MarkerFaceColor','g'); 
               
        end
    end
end
end 
