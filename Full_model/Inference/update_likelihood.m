function [loglklh] = update_likelihood(trials,variational_samples,parameter_current, ...
    background_rate,lklh_func,spike_curves,neurons,prior_info,inference_params,pre_density)
%%
figure_flag =false;
% tic;
% tstart=toc;
% t1=toc;
% ts(:)=0;
n_shape=inference_params.MCsamples_for_gradient;
GP_params=prior_info.prior_parameters.GP_params;
n_cell=length(variational_samples);
n_trial=length(trials);
stim_size=zeros(n_trial,n_cell);
boundary_params = prior_info.prior_parameters.boundary_params;
% 
for  i_trial = 1:n_trial
    for i_loc = 1:size(trials(i_trial).locations,1)
         cell_and_pos=trials(i_trial).cell_and_pos{i_loc};
         if ~isempty(cell_and_pos)   
         power_tmp = trials(i_trial).power_levels(i_loc);
             
           for i=1:size(cell_and_pos,1)
               i_cell = cell_and_pos(i,1);i_pos= cell_and_pos(i,2);
                stim_size(i_trial,i_cell)=stim_size(i_trial,i_cell)+...
                    power_tmp*variational_samples(i_cell).shapes(i_pos);
           end
         end
    end
end
% t2=toc;

Tmax=spike_curves.time_max;Tmin=0;
if isfield(inference_params, 'event_range')
    Tmax = inference_params.event_range(2);
    Tmin = inference_params.event_range(1);
end
max_grid =length(pre_density.pdf_grid);
minimum_stim_threshold=spike_curves.minimum_stim_threshold;
% t3=toc;
loglklh=zeros(n_cell,1);
% Use evenly-spaced grid for spike_curves.current for easy mapping:
current_lb=min(spike_curves.current);
current_gap=spike_curves.current(2)-spike_curves.current(1);
current_max_grid = length(spike_curves.current);

spike_curves_var=spike_curves.sd.^2;
spike_curves_mean=spike_curves.mean;
pdf_grid=pre_density.pdf_grid;
cdf_grid=pre_density.cdf_grid;
grid_bound=pre_density.grid.bound;
grid_gap=pre_density.grid.gap;

gain_sample=reshape([variational_samples(:).gain], [n_cell 1]);
PR_sample=reshape([variational_samples(:).PR], [n_cell 1]);
if isfield(variational_samples(1),'delay_mu')
    delay_mu_sample=reshape([variational_samples(:).delay_mu], [n_cell 1]);
    delay_var_sample=reshape([variational_samples(:).delay_sigma], [n_cell 1]).^2;
else
    delay_mu_sample=zeros([n_cell 1]);
    delay_var_sample=0.001*ones([n_cell 1]).^2;
end
%
% t4=toc;
loglklh_vec = zeros(n_trial,1);
for  i_trial = 1:n_trial
%     t5=toc;
    event_times=trials(i_trial).event_times;
    if  isnan(event_times)
        event_times=[];
    end
%     prob_this_trial=zeros(1,length(event_times)+1);
    prob_this_trial=background_rate*ones(1,length(event_times)+1);
    prob_this_trial(1,end)=background_rate*(Tmax-Tmin);
     
%     t6=toc;
   
    i_count = 1;
    for i_cell = 1:n_cell % can reduce to cell with sufficiently large stimuliaton
        effective_stim=stim_size(i_trial,i_cell)*gain_sample(i_cell);
        if effective_stim>minimum_stim_threshold
            t8=toc;
             delay_mu_temp=delay_mu_sample(i_cell);
        delay_var_temp=delay_var_sample(i_cell);
        i_count = i_count + 1;
        stim_index= min(current_max_grid,...
            max(1,round((effective_stim-current_lb)/current_gap)));
        spike_times_cond_shape=spike_curves_mean(stim_index);
        expectation=delay_mu_temp+spike_times_cond_shape;
        standard_dev=sqrt(delay_var_temp+  mean(spike_curves_var(stim_index)));
%         t9=toc;
        cdf_index_max = max(1,min(max_grid,round( ((Tmax-expectation)/standard_dev +grid_bound)/grid_gap)));
        cdf_index_min = max(1,min(max_grid,round( ((Tmin-expectation)/standard_dev +grid_bound)/grid_gap)));
        pdf_index = max(1,min(max_grid,round( ((event_times-expectation)/standard_dev +grid_bound)/grid_gap)));
        prob_this_trial(i_count,:)=...
            PR_sample(i_cell)*[pdf_grid(pdf_index) cdf_grid(cdf_index_max)-cdf_grid(cdf_index_min)];
%         t10=toc;
%         ts(6)=t9-t8+ts(6);ts(7)=t10-t9+ts(7);
        end
        
    end
%     t11=toc;
    loglklh_vec(i_trial)=  lklh_func(trials(i_trial),prob_this_trial);
%     t7=toc;
%     ts(4)=t6-t5+ts(4);ts(5)=t7-t11+ts(5);
    %%
%     if figure_flag
%         prob_this_trial=zeros(1,Tmax);
%         prob_this_trial(1,:)=background_rate*ones(1,Tmax);
%         for i_cell = 1:n_cell
%             delay_mu_temp=delay_mu_sample(i_cell);
%             delay_sigma_temp=delay_sigma_sample(i_cell);
%             stim_index = ones(n_shape,1);
%             for i_shape = 1:n_shape
%                 stim_temp =stim_size(i_trial,i_cell,i_shape);
%                 effective_stim=stim_temp*gain_sample(i_cell);
%                 if effective_stim>minimum_stim_threshold
%                     stim_index(i_shape)= min(current_max_grid,...
%                         max(1,round((effective_stim-current_lb)/current_gap)));
%                 end
%             end
%             spike_times_cond_shape=spike_curves.mean(stim_index);
%             expectation=delay_mu_temp+mean(spike_times_cond_shape); % mean expectation
%             standard_dev=sqrt(delay_sigma_temp^2+...
%                 mean(spike_curves.sd(stim_index).^2)+ var(spike_times_cond_shape));
%             pdf_index = normpdf(1:Tmax,expectation, standard_dev);
%             prob_this_trial(i_cell,:)=pdf_index;
%         end
%         figure(i_trial)
%         total_prob=sum(prob_this_trial);
%         plot(total_prob)
%         hold on;
%         scatter(event_times, max(total_prob) )
%     end
    
    %%
end
loglklh=sum(loglklh_vec);

% ts(1)=t2-t1+ts(1);ts(2)=t3-t2+ts(2);ts(3)=t4-t3+ts(3);
% tend=toc;


%% Code for integral
% stim_size=zeros(n_trial,n_cell,n_shape);
% for i_cell = 1:n_cell
%     if length(neurons(i_cell).location)==1
%         shifts=[variational_samples(i_cell).shift_x];
%     else
%         shifts=[variational_samples(i_cell).shift_x variational_samples(i_cell).shift_y variational_samples(i_cell).shift_z];
%     end
%     adjusted_location(i_cell,:)=neurons(i_cell).location-shifts;
%     relevant_trials=[];
%     relevant_stim_locations=zeros(0,length(neurons(i_cell).location));
%     trial_and_stim=zeros(0,2);
%     
%     this_adjusted_loc=adjusted_location(i_cell,:);
%     for i_trial = 1:length(trials)
%         relevant_flag = false;
%         for i_loc = 1:size(trials(i_trial).locations,1)
%             rel_position=trials(i_trial).locations(i_loc,:)-this_adjusted_loc;
%             if length(rel_position)==1
%                 this_relevant_flag = (abs(rel_position)<boundary_params);
%             else
%                 this_relevant_flag =check_in_boundary(rel_position,boundary_params);
%             end
%             if  this_relevant_flag
%                 relevant_stim_locations=[relevant_stim_locations; rel_position];
%                 trial_and_stim=[trial_and_stim;  [i_trial i_loc]];
%                 relevant_flag =true;
%             end
%         end
%         if relevant_flag
%             relevant_trials = [relevant_trials i_trial];
%         end
%     end
%     if ~isempty(relevant_trials)
%         if size(relevant_stim_locations,2)==1
%             [GP_samples] = draw_1D_GP(relevant_stim_locations,n_shape,GP_params);
%         else
%             [GP_samples] = draw_3D_GP(relevant_stim_locations,n_shape,GP_params);
%         end
%         for i_shape = 1:n_shape
%             shape_val=GP_samples.full.samples(:,i_shape);
%             for i_rel = 1:length(shape_val)
%                 i_trial = trial_and_stim(i_rel,1);
%                 i_loc = trial_and_stim(i_rel,2);
%                 if length(trials(i_trial).power_levels)==1
%                     power_tmp =  trials(i_trial).power_levels;
%                 else
%                     power_tmp = trials(i_trial).power_levels(i_loc);
%                 end
%                 stim_size(i_trial,i_cell,i_shape)=stim_size(i_trial,i_cell,i_shape)+...
%                     power_tmp*shape_val(i_rel);
%             end
%         end
%         
%         % debugging:
%         %%
%         if figure_flag & (i_trial < 20)
%             figure(i_cell)
%             for i_shape = 1:n_shape
%                 shape_val=GP_samples.full.samples(:,i_shape);
%                 scatter(relevant_stim_locations,shape_val)
%                 hold on;
%             end
%         end
%         %%
%     end
% end


% stim_index = ones(n_shape,1);
%         for i_shape = 1:n_shape
%             stim_temp =stim_size(i_trial,i_cell,i_shape);
%             effective_stim=stim_temp*gain_sample(i_cell);
%             if effective_stim>minimum_stim_threshold
%                 stim_index(i_shape)= min(current_max_grid,...
%                     max(1,round((effective_stim-current_lb)/current_gap)));
%             end
%         end


        
%         spike_times_cond_shape=spike_curves.mean(stim_index);
%         expectation=delay_mu_temp+mean(spike_times_cond_shape); % mean expectation
%         standard_dev=sqrt(delay_sigma_temp^2+...
%             mean(spike_curves.sd(stim_index).^2)+ var(spike_times_cond_shape));
%         cdf_index_max = max(1,min(length(cdf_grid),round( ((Tmax-expectation)/standard_dev +grid.bound)/grid.gap)));
%         cdf_index_min = max(1,min(length(cdf_grid),round( ((Tmin-expectation)/standard_dev +grid.bound)/grid.gap)));
%         pdf_index = max(1,min(length(pdf_grid),round( ((event_times-expectation)/standard_dev +grid.bound)/grid.gap)));
%         prob_this_trial(i_cell,:)=...
%             PR_sample(i_cell)*[pdf_grid(pdf_index) cdf_grid(cdf_index_max)-cdf_grid(cdf_index_min)];
