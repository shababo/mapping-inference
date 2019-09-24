function [parameter_history, elbo_rec] = fit_VI(...
    trials,neurons, ...
    variational_params,prior_params,...
    inference_params,prior_info)
%%
% Update these parameters in parallel 
% par_fields = {'PR' 'delay_var'};
% if ~isfield(variational_params(1),'PR')
%     par_fields(ismember(par_fields,'PR')) = [];
% end
% if ~isfield(variational_params(1),'delay_var')
%     par_fields(ismember(par_fields,'delay_var')) = [];
% end

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
boundary_params = prior_info.GP_params.GP_boundary;
GP_params=prior_info.GP_params;
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
grid.bound=15;grid.gap=0.01;
pre_density.grid=grid;
pre_density.normal_grid = -grid.bound:grid.gap:grid.bound;
for i = 1:length(pre_density.normal_grid)
    pre_density.cdf_grid(i) =normcdf(pre_density.normal_grid(i),0,1);
    pre_density.pdf_grid(i)=normpdf(pre_density.normal_grid(i),0,1);
end
% S=200;
% clear('parameter_history')
% clear('gradients')
%%
tic;
fprintf('Threshold %d; S %d\n',eta_threshold,S)
% ts=zeros(10,1);
while (change_history(iteration) > epsilon && iteration<maxit)
    %%
    tstart=toc;
    parameter_history(iteration,:)=parameter_current;
    % Save the output (parameter_history) to a temp file  
    if mod(iteration, 10)==0
    save('./vi_path_tmp.mat', 'parameter_history');
    end
    iteration=iteration+1;
    adjusted_location=zeros(n_cell,3);
    loglklh=zeros(S,1);logprior=zeros(S,1);logvariational=zeros(S,1);
  
    vsam=cell([S 1]);rsam=cell([S 1]);
   
    for s=1:S
            %         t1=toc;
            [variational_samples,raw_samples] = draw_samples_from_var_dist(parameter_current);
               %---------------------------------%
        % Manually set the true paramters for debugging only:
%         variational_samples(1).gain=0.04;
%         variational_samples.delay_mean=neurons.truth.delay_mean;
%         variational_samples.delay_var=neurons.truth.delay_var;
%         variational_samples.z=1;
%         variational_samples.xy=1;
%         variational_samples(1).background=1e-4;
%   variational_samples(1).PR=rand(1);
%  variational_samples(2).PR=0.4;
% variational_samples.shapes(1)=neurons(1).truth.shapes_sim.values(1);
% variational_samples.shapes(2)=neurons(1).truth.shapes_sim.values(2);
% variational_samples.shapes(2)=rand(1)*(0.9+0.4)-0.4;
        %----------------------------------%
            vsam{s}=variational_samples;rsam{s}=raw_samples;

            logprior(s)=get_logdistribution(variational_samples,raw_samples,prior_params);
            logvariational(s)=get_logdistribution(variational_samples,raw_samples,parameter_current);
            if isfield(variational_samples(1),'background')
                background_rate =  variational_samples(1).background;
            else
                background_rate = eps;
            end
            
            % Calculate likelihood 
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
%         end
    end
%% 
    shapes_val=zeros(S,length(vsam{1}.shapes));
      lklhweights=zeros(S,1);
    for i = 1:S
        shapes_val(i,:) = vsam{i}(1).shapes+vsam{i}(1).gain;
        lklhweights(i)=logprior(i)+loglklh(i)-logvariational(i);
    end
%     for i = 1:size(shapes_val,2)
i=2;
    figure(i)
    scatter(shapes_val(:,i),loglklh,'MarkerFaceColor','blue')
    hold on;
    line(log([neurons.truth.shapes_sim.values(i) neurons.truth.shapes_sim.values(i)]*0.04), [min(loglklh) max(loglklh)],'Color','red','LineWidth',3)
    hold off;
    xlabel('log-Excitability (shape + gain)','FontSize',15);
    ylabel('Log-Likelihood','FontSize',15)
    title('Global','FontSize',15)
%     end
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
