function [] = visualize_fitted_trials(trials, neurons,plot_params)
ref_size=plot_params.ref_size;

if ~isfield(plot_params, 'by_neuron')
    plot_params.by_neuron=false;
end
if plot_params.prediction
    chosen_field = 'predicted';
else
      chosen_field = 'fitted';
end
fig_num = randi(10000);
i = 0;
% figure(fig_num)
tmp_loc=reshape([trials(:).locations],3,[])';
[locations_unique, ~, ic] = unique(tmp_loc,'rows');
% for i_cell = 1:length(neurons)
%     tmp=neurons(i_cell).location;
%     scatter3(tmp(1),tmp(2),tmp(3),ref_size*2,  'MarkerFaceColor','red','MarkerEdgeColor','red');
%     hold on;
%     if plot_params.by_neuron
%         text(tmp(1),tmp(2), tmp(3),num2str(i_cell));
%     end
% end
% 
% for i = 1:size(locations_unique,1)
%     tmp=locations_unique(i,:);
%     % Size scaled with power
%     tmp_size = ref_size/2;
%     scatter3(tmp(1),tmp(2), tmp(3), tmp_size,'MarkerFaceColor','blue','MarkerEdgeColor','blue','MarkerFaceAlpha',0.1);
%     hold on;
%     if ~plot_params.by_neuron
%         text(tmp(1),tmp(2), tmp(3),num2str(i));
%     end
% end
% hold off;

% Find the maximum intensity and use it to standardize all intensities 
plot_params.max_intensity=1e-2;
for i_trial = 1:length(trials)
    this_trial =trials(i_trial);
if plot_params.spike 
    fits = this_trial.(chosen_field).intensity.spike;
else
    fits =this_trial.(chosen_field).intensity.event;
end
    plot_params.max_intensity = max(plot_params.max_intensity, max(max(fits)));
end

% Visualize the trials
if plot_params.by_neuron
    for i_neuron = 1:length(neurons)
       fig= figure(i_neuron+1+fig_num)
        these_trials =trials;
        plot_params.vertical_gap = 0.1;
        stim_size = zeros(length(these_trials),1);
        for i = 1:length(these_trials)
            stim_size(i)=these_trials(i).(chosen_field).stim(i_neuron+1);
        end
%         [~, tmp]=unique(stim_size); % Ascending
%         
%         rankings = 1:length(these_trials);
%         rankings(tmp)=rankings;
        visualize_fitted_trial_multiple(these_trials, stim_size,plot_params)
        title(['Neuron ' num2str( i_neuron )])
        if isfield(plot_params,'save_path')
            save_path = [plot_params.save_path '_' num2str(i+1) '.png'];
            saveas(fig,save_path)
        end
    end
else
    for i = 1:size(locations_unique,1)

        these_trials =trials(ic==i);

                
                
%         fig=figure(i+1+fig_num);
% %         plot_params.vertical_gap = 0.1;
% % %         [~, tmp]=sort([these_trials(:).power_levels]); % Ascending
% % %         rankings = 1:length(these_trials);
% % %         rankings(tmp)=rankings;        
%         visualize_fitted_trial_multiple(these_trials,[these_trials(:).power_levels],plot_params)
%         title(['Location ' num2str(i)])
%         if isfield(plot_params,'save_path')
%             save_path = [plot_params.save_path '_' num2str(i+1) '.png'];
%             saveas(fig,save_path)
%         end
        
        the_trial = find([these_trials.power_levels] == 75);
        the_trial = these_trials(the_trial(1))
        fig=figure(199999);%figure(size(locations_unique,1)+fig_num+200);
        hold on
        scatter3(locations_unique(i,1),locations_unique(i,2),locations_unique(i,3),[],the_trial.fitted.stim(2)/75)
        
        fig=figure(1999999);%figure(size(locations_unique,1)+fig_num+201);
        if isnan(the_trial.event_times) || the_trial.event_times > 160
            event_time = 0;
        else
            event_time = (160 - the_trial.event_times)/160;
        end
        hold on
        scatter3(locations_unique(i,1),locations_unique(i,2),locations_unique(i,3),[],event_time)

        
    end
end



