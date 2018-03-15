function [loglklh] = update_likelihood_dev3(...
    gamma_sample_mat,gain_sample_mat,stim_size,mpp,...
    delay_mu_sample_mat,delay_sigma_sample_mat,...
    prob_trace_full,minimum_stim_threshold,stim_scale,background_rate,...
    relevant_trials,lklh_func,spike_curves)
%%
n_cell=size(gamma_sample_mat,1);
S=size(gamma_sample_mat,2);
n_trial=length(mpp);

Tmax=300; % need to fix this 
% Calculate a grid for standard normal r.v. for quick CDF calculation:
grid.bound=4;
grid.gap=0.1;
normal_grid = -grid.bound:grid.gap:grid.bound;
for i = 1:length(normal_grid)
    cdf_grid(i) = normcdf(normal_grid(i),0,1);
end
loglklh=zeros(n_cell,S);
% tic;
% time_rec=zeros(10,1);

% Check the likelihood?

for s=1:S
    gamma_sample=gamma_sample_mat(:,s);gain_sample=gain_sample_mat(:,s);
    delay_mu_sample=delay_mu_sample_mat(:,s);delay_sigma_sample=delay_sigma_sample_mat(:,s);
%     delay_mu_sample=0;delay_sigma_sample=0;
%%
    loglklh_vec = zeros(n_trial,1);
    for  i_trial = 1:n_trial
%          ts=toc;

event_times=mpp(i_trial).event_times;
stim_temp =stim_size(i_trial,:)';
effective_stim= stim_temp.*gain_sample;
stimulated_cells = find(effective_stim>minimum_stim_threshold);
effective_stim=effective_stim(stimulated_cells );
delay_mu_temp=delay_mu_sample(stimulated_cells);
delay_sigma_temp=delay_sigma_sample(stimulated_cells);
%         te=toc;time_rec(1)=time_rec(1)+te-ts;

% ts=toc;
% find the index of the actual stimulation in the vector:
if ~isempty(stimulated_cells)
    stim_index=zeros(length(stimulated_cells),1);
    for i_stim = 1:length(stimulated_cells)
        [~, stim_index(i_stim)]=min(abs(effective_stim(i_stim) - spike_curves.current));
    end
end

%              te=toc;time_rec(2)=time_rec(2)+te-ts;

%stim_index=min(n_grid,max(1,round(effective_stim*stim_scale)));
%       prob_this_trial= (gamma_sample(stimulated_cells)*ones(1,n_grid)).*prob_trace_full(stim_index,:);
%       prob_this_trial=[background_rate*ones(1, size(prob_this_trial,2)); prob_this_trial];

if ~isempty(stimulated_cells)
    prob_this_trial=zeros(length(stimulated_cells),length(event_times)+1);
    % the extra one is for no event
    
    for i_stim = 1:length(stimulated_cells)
        
        
        expectation=delay_mu_temp(i_stim)+spike_curves.mean(stim_index(i_stim));
        standard_dev=sqrt(delay_sigma_temp(i_stim)^2+...
            spike_curves.sd(stim_index(i_stim))^2+spike_curves.dev(stim_index(i_stim))^2);
        %                         ts=toc;
        cdf_index = max(1,min(length(cdf_grid),round( ((Tmax-expectation)/standard_dev +grid.bound)/grid.gap)));
        prob_this_trial(i_stim,:)=gamma_sample(i_stim)*...
            [normpdf(event_times,expectation,standard_dev) cdf_grid(cdf_index)];
        %                     normcdf(Tmax,expectation,standard_dev)];
        %             te=toc;
        % time_rec(3)=time_rec(3)+te-ts;
        % prob_test=normpdf(1:Tmax,expectation,standard_dev);
    end
    
    prob_this_trial(end,:)=background_rate*ones(1,length(event_times)+1);
    prob_this_trial(end,end)=background_rate*Tmax; % this is an approximation
    
else
    prob_this_trial=background_rate*ones(1,length(event_times)+1);
    prob_this_trial(end)=background_rate*Tmax;
end

%         plot(1:300,prob_this_trial(1,:))
%         hold on;
%         scatter(mpp(i_trial).event_times,0)
        
% ts=toc;
        [loglklh_vec(i_trial)]=  lklh_func(mpp(i_trial),prob_this_trial);
%             te=toc;time_rec(4)=time_rec(4)+te-ts;

    end
%     ts=toc;
    for i_cell = 1:n_cell
        loglklh(i_cell,s)=sum(loglklh_vec(relevant_trials{i_cell}));
    end
%          te=toc;time_rec(5)=time_rec(5)+te-ts;
%%
end

%%
% scatter(gain_sample_mat,loglklh)
