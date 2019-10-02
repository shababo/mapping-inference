function [] = visualize_fitted_trials(trials, neurons,plot_params)
%%
ref_size=plot_params.ref_size;

if ~isfield(plot_params, 'by_neuron')
    plot_params.by_neuron=false;
end
if plot_params.prediction
    chosen_field = 'predicted';
else
      chosen_field = 'fitted';
end
if ~isfield(plot_params, 'stim_scale')
     plot_params.stim_scale = @(x) x;
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
        plot_params.gap = 0.1;
        plot_params.colors=lines(size(locations_unique,1)); 
        plot_params.loc_indices=ic;
        stim_size = zeros(length(these_trials),1);
        
        for i = 1:length(these_trials)
            stim_size(i)=these_trials(i).(chosen_field).stim(i_neuron+1);
        end
        stim_size = log(stim_size);
%         [~, tmp]=unique(stim_size); % Ascending
%         
%         rankings = 1:length(these_trials);
%         rankings(tmp)=rankings;

        plot_params.these_indices=plot_params.loc_indices;
        [covered_flags]=visualize_fitted_trial_multiple(these_trials, stim_size,plot_params);
        title(plot_params.main_title,'FontSize', plot_params.lab_size)
        if isfield(plot_params,'save_path')
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
            save_path = [plot_params.save_path '_All.png'];
            saveas(fig,save_path)
        end
    end
else
    for i = 1:size(locations_unique,1)
these_trials =trials(ic==i);
        fig=figure(i+1+fig_num);
%         plot_params.vertical_gap = 0.1;
% %         [~, tmp]=sort([these_trials(:).power_levels]); % Ascending
% %         rankings = 1:length(these_trials);
% %         rankings(tmp)=rankings;        

        plot_params.these_indices=plot_params.loc_indices(ic==i);
        visualize_fitted_trial_multiple(these_trials,[these_trials(:).power_levels],plot_params)
        title(['Location ' num2str(i)],'FontSize', plot_params.lab_size)
        if isfield(plot_params,'save_path')
            save_path = [plot_params.save_path '_' num2str(i+1) '.png'];
            saveas(fig,save_path)
        end
    end
end

% set(gca, 'XScale',  plot_params.stim_scale)
% xlim([-3 5])
yyaxis left
ylim([0 plot_params.time_max]);

yyaxis right
ylabel('Expected # events','FontSize', plot_params.lab_size)
ylim([0 plot_params.time_max]);

% if isfield(plot_params,'off_loc')
%     fig= figure(i_neuron+4+fig_num)
%        these_trials =trials;
%        plot_params.vertical_gap = 0.1;
%        dist_list = zeros(length(these_trials),1);
%         for i = 1:length(these_trials)
%             
%             rel_pos=these_trials(i).locations-neurons(1).location;
%             dist_list(i)= min( sqrt(sum(rel_pos(1:2).^2)),abs(rel_pos(3)));
%         end
%        visualize_fitted_trial_multiple(these_trials,dist_list,plot_params)
%       ylim([0,2.5])
%      ylabel('Trials ordered by deviation from xy-plane and z-axis');
% end
% 

%% Draw a spatial map: 
fig=figure(10+fig_num);
for i = 1:size(locations_unique,1)
    scatter3(locations_unique(i,1),locations_unique(i,2),locations_unique(i,3),plot_params.markerSize,'filled',...
        'MarkerFaceColor', plot_params.colors(i,:),  'MarkerEdgeColor',plot_params.colors(i,:))
    hold on;
    if ~covered_flags(i)
       text(locations_unique(i,1),locations_unique(i,2),locations_unique(i,3), ['Loc ' num2str(i)]);
    end
    if isfield(plot_params,'save_path')
        save_path = [plot_params.save_path '_Map.png'];
        saveas(fig,save_path)
    end
end
%% Draw the trials at the problematic locations: 

plot_params.fit_type='full_intensity';
for i = 1:size(locations_unique,1)
    if ~covered_flags(i)
        these_trials =trials(ic==i);
        fig=figure(10+i+fig_num);
        plot_params.these_indices=plot_params.loc_indices(ic==i);
        %         plot_params.vertical_gap = 0.1;
        % %         [~, tmp]=sort([these_trials(:).power_levels]); % Ascending
        % %         rankings = 1:length(these_trials);
        % %         rankings(tmp)=rankings;
        visualize_fitted_trial_multiple(these_trials,[these_trials(:).power_levels],plot_params);
        title(['Location ' num2str(i)],'FontSize', plot_params.lab_size)
            if isfield(plot_params,'save_path')
                save_path = [plot_params.save_path '_Loc' num2str(i) '.png'];
                saveas(fig,save_path)
            end
        ylim([0 5+plot_params.time_max]);
    end
end
